import streamlit as st
import subprocess
import os
import re
import time
import threading
import pandas as pd

# ── ANSI / log helpers ────────────────────────────────────────────────────────
_ANSI_RE = re.compile(r'\033\[[0-9;]*m')
_OPT_TYPES = {"START", "PROGRESS", "RESULT", "WARNING", "ERROR", "DONE"}


def _strip_ansi(s: str) -> str:
    return _ANSI_RE.sub('', s)


def _fmt_line(raw: str) -> str:
    """Parse one opt log line and return an HTML-formatted string with colors."""
    clean = _strip_ansi(raw).strip()
    if not clean:
        return ""
    parts = clean.split(None, 3)
    if len(parts) == 4 and parts[2] in _OPT_TYPES:
        mtype, msg = parts[2], parts[3]
        if mtype == "START":
            return f'<span style="color:#1f77b4;font-weight:bold">&#9654; {msg}</span>'
        elif mtype == "PROGRESS":
            return f'<span style="color:#888888">&mdash; {msg}</span>'
        elif mtype == "RESULT":
            return f'<span style="color:#2ca02c;font-weight:bold">&rarr; {msg}</span>'
        elif mtype == "WARNING":
            return f'<span style="color:#e07b00;font-weight:bold">&#9888; {msg}</span>'
        elif mtype == "ERROR":
            return f'<span style="color:#d62728;font-weight:bold">&#10007; Error: {msg}</span>'
        elif mtype == "DONE":
            return f'<span style="color:#2ca02c">&#10003; {msg}</span>'
    return f'<code style="color:#555555;font-size:0.85em">{clean}</code>'


def browse_file(session_key: str, prompt: str = "Select a file") -> None:
    """Open a native macOS file-picker via osascript (AppleScript).

    Avoids all tkinter/thread issues — osascript manages its own event loop.
    session_key must NOT be the key= of any widget; use value= instead.
    """
    script = f'POSIX path of (choose file with prompt "{prompt}")'
    try:
        result = subprocess.run(
            ["osascript", "-e", script],
            capture_output=True, text=True, timeout=120,
        )
        path = result.stdout.strip()
        if path:
            st.session_state[session_key] = path
            st.rerun()
    except subprocess.TimeoutExpired:
        st.warning("File browser timed out. Please type the path manually.")
    except Exception as e:
        st.warning(f"Could not open file browser: {e}. Please type the path manually.")


# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="OPT — Off-target Probe Tracker",
    page_icon="🧬",
    layout="wide",
)

DEFAULT_SCHEMA = ["transcript", "ID", "Parent", "gene_name", "transcript_type"]

SCHEMA_PRESETS = {
    "GENCODE": "transcript,ID,Parent,gene_name,transcript_type",
    "RefSeq":  "transcript,ID,Parent,gene,gbkey",
    "CHESS":   "transcript,ID,Parent,gene_name,gene_type",
}

_REF_COLORS = {"GENCODE": "#1f77b4", "CHESS": "#17becf", "RefSeq": "#2ca02c"}
_MULTI_PRESET = "All (GENCODE + CHESS + RefSeq)"


# ── Biotype color helpers ─────────────────────────────────────────────────────
def _normalize_biotype(biotype: str) -> str:
    """Return canonical group name used for filtering and sorting."""
    if biotype in ("protein_coding", "mRNA"):
        return "protein coding"
    if "pseudogene" in biotype:
        return "pseudogene"
    if biotype in ("lncRNA", "ncRNA"):
        return "lncRNA / ncRNA"
    return biotype.replace("_", " ")


def _biotype_color(biotype: str) -> str:
    if biotype in ("protein_coding", "mRNA"):
        return "#d62728"
    if "pseudogene" in biotype:
        return "#888888"
    if biotype in ("lncRNA", "ncRNA"):
        return "#e07b00"
    if biotype in ("miscRNA", "miRNA", "snoRNA", "snRNA"):
        return "#9467bd"
    return "#1f77b4"


def _biotype_badge(biotype: str) -> str:
    color = _biotype_color(biotype)
    # Normalize label for grouped types; pseudogenes keep their specific name
    if biotype in ("protein_coding", "mRNA"):
        label = "protein coding"
    elif biotype in ("lncRNA", "ncRNA"):
        label = "lncRNA / ncRNA"
    else:
        label = biotype.replace("_", " ")
    return (
        f'<span style="background:{color};color:white;padding:2px 7px;'
        f'border-radius:3px;font-size:0.78em;margin:1px;display:inline-block">'
        f'{label}</span>'
    )


# ── Session state ─────────────────────────────────────────────────────────────
def init_session_state():
    defaults = {
        "run_stdout":     None,
        "run_returncode": None,
        "run_module":     None,
        "out_dir":        None,
        "last_cmd":       None,
        "current_step":   0,
        "biotype_filter": "All",
        "table_sort_col": "Target gene",
        "table_sort_asc": True,
        "opt_running":         False,
        "opt_log_path":        None,
        "opt_rc_path":         None,
        "opt_multi_mode":      False,
        "opt_runs_queue":      [],
        "opt_current_run":     {},
        "opt_completed_runs":  [],
        "opt_out_dir_base":       None,
        "run_annotation_preset":  "",
        # Path fields managed via Browse (not owned by any widget via key=)
        "all_query_path_val":     "",
        "all_target_val":         "",
        "all_annotation_val":     "",
        "all_syn_val":            "",
        "all_out_dir_val":        "./opt_results",
        "multi_gencode_ann_val":  "",
        "multi_gencode_fa_val":   "",
        "multi_chess_ann_val":    "",
        "multi_chess_fa_val":     "",
        "multi_refseq_ann_val":   "",
        "multi_refseq_fa_val":    "",
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


# ── Header ────────────────────────────────────────────────────────────────────
def render_header() -> None:
    st.title("OPT — Off-target Probe Tracker")
    st.markdown(
        """
        OPT identifies potential off-target binding of probe sequences against a reference transcriptome using nucleotide alignment (**nucmer**). The goal of OPT is to help evaluate probe specificity before experiments by detecting probes that may hybridize to unintended transcripts.

        The pipeline runs three steps automatically:

        **flip** - aligns probe sequences to the source transcriptome to determine the correct strand orientation and ensures probes are evaluated in the proper direction.

        **track** - aligns the corrected probe sequences to all transcripts in the selected reference transcriptome and identifies probes that map to more than one gene, flagging these as potential off-target probes.

        **stat** - aggregates the alignment results and summarizes off-target binding events by target gene, off-target gene, and transcript biotype.

        ### Using this app

        This Streamlit application provides an interactive interface to run OPT without using the command line. Users can:

        -  Load probe sequences and reference annotation files  
        -  Configure alignment parameters  
        -  Run the OPT pipeline directly from the browser  

        ### Output and results

        After the analysis completes, the app summarizes results with:

        -  A table listing each **target gene**, the number of **probes with detected off-targets**, the corresponding **off-target genes**, 
        the **gene biotype** of the off-target transcripts, the **CIGAR alignment patterns** describing how probes aligned to those transcripts,
        and the reference annotation.

        - View probe-level predicted off-target details to inspect individual probes and their alignments.

        The Streamlit interface highlights the key off-target results so users can quickly identify problematic probes. All original OPT output files are also saved to the specified output directory for further inspection or downstream analysis.

        If you find this tool useful in your research, please consider citing:

        Hallinan et al., *eLife* 2025  
        https://elifesciences.org/reviewed-preprints/107070
        """
    )


# ── Reusable path field (text input + Browse + real-time validation) ───────────
def _path_field(
    caption: str,
    state_key: str,
    browse_key: str,
    browse_prompt: str,
    placeholder: str,
    help_text: str = "",
    optional: bool = False,
) -> str:
    label = caption + (" (optional)" if optional else " *")
    st.caption(label)
    col_in, col_btn = st.columns([5, 1])
    with col_in:
        path = st.text_input(
            caption,
            value=st.session_state.get(state_key, ""),
            placeholder=placeholder,
            help=help_text,
            label_visibility="collapsed",
        )
        st.session_state[state_key] = path
    with col_btn:
        if st.button("Browse", key=browse_key, use_container_width=True):
            browse_file(state_key, prompt=browse_prompt)
    # Real-time validation feedback
    if path.strip():
        if os.path.exists(path.strip()):
            st.markdown(
                '<span style="color:#2ca02c;font-size:0.82em">&#10003; File found</span>',
                unsafe_allow_html=True,
            )
        else:
            st.markdown(
                '<span style="color:#d62728;font-size:0.82em">&#10007; File not found</span>',
                unsafe_allow_html=True,
            )
    # Show help text as an info box for optional fields (label is collapsed so tooltip is hidden)
    if optional and help_text:
        st.info(help_text)
    return path


# ── Inputs ────────────────────────────────────────────────────────────────────
def render_all_inputs() -> tuple:
    """Returns (global_args dict, module_inputs dict)."""

    # Output dir + threads — always visible at top
    st.subheader("Run Configuration")
    od_col, th_col = st.columns([3, 2])
    with od_col:
        out_dir = st.text_input(
            "Output Directory *",
            value=st.session_state.get("all_out_dir_val", "./opt_results"),
            help="Where OPT writes all output files.",
        )
        st.session_state["all_out_dir_val"] = out_dir
    with th_col:
        threads = st.slider(
            "Threads", min_value=1, max_value=16, value=1,
            help="Number of CPU threads passed to nucmer (-t). More threads speeds up alignment."
        )

    st.divider()
    st.subheader("Input Files")

    query_path = _path_field(
        "Probe FASTA", "all_query_path_val", "all_query_browse",
        "Select Probe FASTA",
        "/path/to/probes.fa",
        "Path to probe sequences (.fa or .fasta). Header format: >gene_id|gene_name|accession",
    )
    st.info(
        "**Probe FASTA header format required:** `>gene_id|gene_name|accession`  \n"
        "Example: `>ENSG00000170458|CD14|22f9405`"
    )

    # Annotation format selection — determines single vs. multi-annotation inputs
    st.caption("Annotation format")
    schema_preset_early = st.radio(
        "Annotation format",
        ["GENCODE", "RefSeq", "CHESS", "Other", _MULTI_PRESET],
        horizontal=True,
        label_visibility="collapsed",
        key="schema_preset_radio",
        help=(
            "Select the source of your annotation GFF file. "
            "Choose 'All' to run against GENCODE, CHESS, and RefSeq simultaneously."
        ),
    )

    multi_mode = schema_preset_early == _MULTI_PRESET

    if not multi_mode:
        target = _path_field(
            "Target transcript FASTA", "all_target_val", "all_target_browse",
            "Select Target Transcript FASTA",
            "/path/to/transcripts.fa",
            "Local path to target transcript sequences (.fa or .fasta).",
        )
        annotation = _path_field(
            "Annotation GFF/GTF", "all_annotation_val", "all_annotation_browse",
            "Select Annotation GFF/GTF",
            "/path/to/annotation.gff",
            "Local path to transcript annotation (.gff, .gff3, or .gtf).",
        )
    else:
        target = ""
        annotation = ""
        st.caption("Provide transcript FASTA and annotation GFF for each reference:")
        for ref_name, fa_key, ann_key in [
            ("GENCODE", "multi_gencode_fa_val", "multi_gencode_ann_val"),
            ("CHESS",   "multi_chess_fa_val",   "multi_chess_ann_val"),
            ("RefSeq",  "multi_refseq_fa_val",  "multi_refseq_ann_val"),
        ]:
            color = _REF_COLORS[ref_name]
            st.markdown(
                f'<span style="background:{color};color:white;padding:2px 8px;'
                f'border-radius:3px;font-size:0.85em;font-weight:bold">{ref_name}</span>',
                unsafe_allow_html=True,
            )
            _path_field(
                f"{ref_name} transcript FASTA", fa_key, f"{ref_name.lower()}_fa_browse",
                f"Select {ref_name} Transcript FASTA",
                "/path/to/transcripts.fa", "",
            )
            _path_field(
                f"{ref_name} annotation GFF", ann_key, f"{ref_name.lower()}_ann_browse",
                f"Select {ref_name} Annotation GFF",
                "/path/to/annotation.gff", "",
            )

    syn_file = _path_field(
        "Gene synonyms CSV", "all_syn_val", "all_syn_browse",
        "Select Gene Synonyms CSV",
        "/path/to/synonyms.csv",
        (
            "Use this if the gene names in your probe FASTA headers differ from those in the "
            "reference annotation — for example, WARS vs WARS1, or MYCN vs N-Myc. "
            "The file should be a two-column CSV with no header: "
            "column 1 = probe gene name, column 2 = annotation gene name. "
            "Example row: WARS,WARS1"
        ),
        optional=True,
    )

    st.divider()
    st.subheader("Analysis Options")

    col_pl_en, col_pl_val = st.columns([1, 2])
    with col_pl_en:
        use_pl = st.checkbox(
            "Pad length (-pl)", value=False,
            help=(
                "Number of bases at each end of the probe where mismatches are allowed. "
                "This parameter was used in the original paper because Xenium probes are circular "
                "and can tolerate mismatches at the termini, while correct base pairing in the "
                "central region is required for ligation."
            ),
        )
    with col_pl_val:
        pl_val = st.number_input(
            "Pad length value",
            min_value=1, max_value=20, value=10, step=1,
            disabled=not use_pl,
            label_visibility="collapsed",
        )
    pad_length = int(pl_val) if use_pl else 0

    col_mm_en, col_mm_val = st.columns([1, 2])
    with col_mm_en:
        use_mm = st.checkbox(
            "Max mismatches anywhere (-mm)", value=False,
            help=(
                "Allow up to N mismatches anywhere in the full probe sequence. "
                "Can be combined with Pad length: when both are set, both conditions "
                "must be satisfied (NM ≤ N AND mismatches confined to terminal pad bases)."
            ),
        )
    with col_mm_val:
        mm_val = st.number_input(
            "Max mismatches value",
            min_value=0, max_value=10, value=5, step=1,
            disabled=not use_mm,
            label_visibility="collapsed",
        )
    max_mismatches = int(mm_val) if use_mm else -1

    # Schema fields — only for single-annotation mode
    if not multi_mode:
        if schema_preset_early == "Other":
            st.caption("GFF/GTF Schema (--schema)")
            st.caption(
                "Five fields used to parse the annotation file. "
                "See the README for guidance on non-standard annotation formats."
            )
            labels = [
                "Feature type (3rd col)",
                "Transcript ID attribute",
                "Parent attribute",
                "Gene name attribute",
                "Transcript type attribute",
            ]
            schema_fields = [
                st.text_input(f"Field {i+1}: {label}", value=DEFAULT_SCHEMA[i], key=f"schema_{i}")
                for i, label in enumerate(labels)
            ]
        else:
            schema_fields = SCHEMA_PRESETS[schema_preset_early].split(",")
    else:
        schema_fields = DEFAULT_SCHEMA  # unused in multi-mode

    with st.expander("Advanced Options", expanded=False):
        bam = st.checkbox(
            "Store alignments as BAM (--bam)", value=False,
            help="Save nucmer alignments as BAM instead of SAM.",
        )
        min_exact = st.number_input(
            "Min exact match length (-l)", min_value=1, value=20, step=1,
            help="Minimum exact match length used by nucmer. Increasing this value speeds up the search but may miss off-targets that share only short exact matches with the probe. Decreasing it increases sensitivity and may detect more potential off-targets, but will increase runtime and may introduce more false positives. The default is 20, which provides a good balance for probes 25-50 nt in length.",
        )
        st.info(
            "For Bowtie2 or a custom aligner, use the CLI version:  \n"
            "https://github.com/JEFworks-Lab/off-target-probe-tracker"
        )

    global_args = {
        "out_dir":       out_dir,
        "threads":       threads,
        "bam":           bam,
        "min_exact":     min_exact,
        "schema":        ",".join(schema_fields),
        "schema_fields": schema_fields,
    }
    module_inputs = {
        "query":           query_path,
        "target":          target,
        "annotation":      annotation,
        "syn_file":        syn_file,
        "pad_length":      pad_length,
        "max_mismatches":  max_mismatches,
        "multi_mode":        multi_mode,
        "annotation_preset": schema_preset_early,
        "multi_refs": {
            "GENCODE": {
                "annotation": st.session_state.get("multi_gencode_ann_val", ""),
                "target":     st.session_state.get("multi_gencode_fa_val",  ""),
                "schema":     SCHEMA_PRESETS["GENCODE"],
            },
            "CHESS": {
                "annotation": st.session_state.get("multi_chess_ann_val", ""),
                "target":     st.session_state.get("multi_chess_fa_val",  ""),
                "schema":     SCHEMA_PRESETS["CHESS"],
            },
            "RefSeq": {
                "annotation": st.session_state.get("multi_refseq_ann_val", ""),
                "target":     st.session_state.get("multi_refseq_fa_val",  ""),
                "schema":     SCHEMA_PRESETS["RefSeq"],
            },
        },
    }
    return global_args, module_inputs


# ── Build command ─────────────────────────────────────────────────────────────
def build_command(global_args: dict, module_inputs: dict) -> list:
    cmd = ["opt"]
    cmd += ["-o", global_args["out_dir"]]
    cmd += ["-p", str(global_args["threads"])]
    if global_args["bam"]:
        cmd.append("--bam")
    cmd += ["-l", str(global_args["min_exact"])]
    cmd += ["--schema", global_args["schema"]]
    cmd.append("--force")

    cmd.append("all")

    cmd += ["-q", module_inputs["query"].strip()]
    cmd += ["-t", module_inputs["target"].strip()]
    cmd += ["-a", module_inputs["annotation"].strip()]
    cmd += ["-pl", str(module_inputs["pad_length"])]

    if module_inputs.get("max_mismatches", -1) >= 0:
        cmd += ["-mm", str(module_inputs["max_mismatches"])]

    syn = module_inputs.get("syn_file", "").strip()
    if syn:
        cmd += ["-s", syn]

    return cmd


# ── Validate inputs ───────────────────────────────────────────────────────────
def _check_ext(path: str, allowed: tuple, label: str, errors: list) -> None:
    """Append an error if path does not end with one of the allowed extensions."""
    if path and not any(path.lower().endswith(ext) for ext in allowed):
        errors.append(f"{label} must be {' or '.join(allowed)} (got: '{path}')")


def validate_inputs(global_args: dict, module_inputs: dict) -> list:
    errors = []
    if not global_args["out_dir"]:
        errors.append("Output directory is required.")

    _FA_EXTS  = (".fa", ".fasta")
    _ANN_EXTS = (".gff", ".gff3", ".gtf")

    if module_inputs.get("multi_mode"):
        # Multi-annotation: validate probe FASTA + all 6 reference files
        query = module_inputs.get("query", "").strip()
        if not query:
            errors.append("Probe FASTA path is required.")
        else:
            _check_ext(query, _FA_EXTS, "Probe FASTA", errors)
            if not os.path.exists(query):
                errors.append(f"Probe FASTA does not exist: '{query}'")
        for ref_name, ref_info in module_inputs.get("multi_refs", {}).items():
            for ftype, key, exts in [
                ("annotation GFF", "annotation", _ANN_EXTS),
                ("transcript FASTA", "target",   _FA_EXTS),
            ]:
                path = ref_info.get(key, "").strip()
                if not path:
                    errors.append(f"{ref_name} {ftype} path is required.")
                else:
                    _check_ext(path, exts, f"{ref_name} {ftype}", errors)
                    if not os.path.exists(path):
                        errors.append(f"{ref_name} {ftype} does not exist: '{path}'")
    else:
        for field in global_args.get("schema_fields", []):
            if "," in field:
                errors.append(f"Schema field contains a comma: '{field}'. Remove it.")
        for field, label, exts in [
            ("query",      "Probe FASTA",            _FA_EXTS),
            ("target",     "Target transcript FASTA", _FA_EXTS),
            ("annotation", "Annotation GFF/GTF",      _ANN_EXTS),
        ]:
            path = module_inputs.get(field, "").strip()
            if not path:
                errors.append(f"{label} path is required.")
            else:
                _check_ext(path, exts, label, errors)
                if not os.path.exists(path):
                    errors.append(f"{label} path does not exist: '{path}'")

    syn = module_inputs.get("syn_file", "").strip()
    if syn and not os.path.exists(syn):
        errors.append(f"Gene synonyms file does not exist: '{syn}'")
    return errors


# ── Step progress indicator ───────────────────────────────────────────────────
def _step_html(current_step: int) -> str:
    """Return HTML for a flip → track → stat progress bar.
    current_step: 0=idle, 1=flip running, 2=track running, 3=stat running, 4=done.
    """
    steps = ["flip", "track", "stat"]
    parts = []
    for i, step in enumerate(steps):
        n = i + 1
        if n < current_step:  # completed (also handles current_step=4: all done)
            icon  = '<span style="color:#2ca02c">&#10003;</span>'
            label = f'<span style="color:#2ca02c">{step}</span>'
        elif n == current_step:
            icon  = '<span style="color:#1f77b4;font-weight:bold">&#9654;</span>'
            label = f'<span style="color:#1f77b4;font-weight:bold">{step}</span>'
        else:
            icon  = '<span style="color:#aaaaaa">&#9675;</span>'
            label = f'<span style="color:#aaaaaa">{step}</span>'
        parts.append(f'{icon} {label}')
    return ' &nbsp;&rarr;&nbsp; '.join(parts)


# ── Run OPT ───────────────────────────────────────────────────────────────────
def _wait_for_proc(proc, rc_path: str) -> None:
    """Daemon thread: waits for proc to finish, writes exit code to rc_path."""
    proc.wait()
    with open(rc_path, "w") as f:
        f.write(str(proc.returncode))


def _start_run(spec: dict) -> None:
    """Launch one opt run as a non-blocking background process."""
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    log_f = open(spec["log_path"], "w")
    proc  = subprocess.Popen(
        spec["cmd"], stdout=log_f, stderr=subprocess.STDOUT, text=True, env=env,
    )
    threading.Thread(target=_wait_for_proc, args=(proc, spec["rc_path"]),
                     daemon=True).start()
    st.session_state["opt_log_path"] = spec["log_path"]
    st.session_state["opt_rc_path"]  = spec["rc_path"]


def build_multi_commands(global_args: dict, module_inputs: dict) -> list:
    """Return a list of run spec dicts (one per reference annotation)."""
    base_out = global_args["out_dir"]
    specs = []
    for ref_name, ref_info in module_inputs["multi_refs"].items():
        sub_dir = os.path.join(base_out, ref_name.lower())
        cmd = [
            "opt",
            "-o", sub_dir,
            "-p", str(global_args["threads"]),
            "-l", str(global_args["min_exact"]),
            "--schema", ref_info["schema"],
            "--force",
            "all",
            "-q", module_inputs["query"].strip(),
            "-t", ref_info["target"].strip(),
            "-a", ref_info["annotation"].strip(),
            "-pl", str(module_inputs["pad_length"]),
        ]
        if module_inputs.get("max_mismatches", -1) >= 0:
            cmd += ["-mm", str(module_inputs["max_mismatches"])]
        syn = module_inputs.get("syn_file", "").strip()
        if syn:
            cmd += ["-s", syn]
        specs.append({"name": ref_name, "cmd": cmd, "out_dir": sub_dir})
    return specs


def _merge_multi_results(completed_runs: list, out_dir: str) -> None:
    """Merge TSV results from all annotation subdirs into out_dir."""
    probe_dfs    = []
    all_dfs      = []
    summary_dfs  = []
    probe_ids_union: set = set()

    for run in completed_runs:
        sub_dir  = run["out_dir"]
        ref_name = run["name"]
        p_path   = os.path.join(sub_dir, "probe2targets_offtargets.tsv")
        a_path   = os.path.join(sub_dir, "probe2targets.tsv")
        s_path   = os.path.join(sub_dir, "collapsed_summary_offtargets.tsv")
        if os.path.exists(p_path):
            df = pd.read_csv(p_path, sep="\t")
            df["reference_annotation"] = ref_name
            probe_dfs.append(df)
            if "probe_id" in df.columns:
                probe_ids_union.update(df["probe_id"].astype(str))
        if os.path.exists(a_path):
            df = pd.read_csv(a_path, sep="\t")
            df["reference_annotation"] = ref_name
            all_dfs.append(df)
        if os.path.exists(s_path):
            df = pd.read_csv(s_path, sep="\t")
            df["reference_annotation"] = ref_name
            summary_dfs.append(df)

    if probe_dfs:
        pd.concat(probe_dfs, ignore_index=True).to_csv(
            os.path.join(out_dir, "probe2targets_offtargets.tsv"), sep="\t", index=False)
    if all_dfs:
        pd.concat(all_dfs, ignore_index=True).to_csv(
            os.path.join(out_dir, "probe2targets.tsv"), sep="\t", index=False)
    if summary_dfs:
        merged = pd.concat(summary_dfs, ignore_index=True)
        if "target_gene" in merged.columns:
            merged = merged.drop_duplicates(subset=["target_gene"], keep="first")
        merged.to_csv(
            os.path.join(out_dir, "collapsed_summary_offtargets.tsv"), sep="\t", index=False)
    with open(os.path.join(out_dir, "stat_off_target_probes.txt"), "w") as f:
        for pid in sorted(probe_ids_union):
            f.write(pid + "\n")


def run_opt(global_args: dict, module_inputs: dict) -> None:
    errors = validate_inputs(global_args, module_inputs)
    if errors:
        for e in errors:
            st.error(e)
        return

    if module_inputs.get("multi_mode"):
        _run_multi_opt(global_args, module_inputs)
        return

    cmd     = build_command(global_args, module_inputs)
    out_dir = global_args["out_dir"]
    os.makedirs(out_dir, exist_ok=True)

    # Persist annotation preset to disk so it survives page refresh / app restart
    preset = module_inputs.get("annotation_preset", "")
    if preset in _REF_COLORS:
        with open(os.path.join(out_dir, "_opt_annotation_preset.txt"), "w") as _f:
            _f.write(preset)

    log_path = os.path.join(out_dir, "_opt_run.log")
    rc_path  = os.path.join(out_dir, "_opt_run.rc")
    for p in [log_path, rc_path]:
        if os.path.exists(p):
            os.remove(p)

    try:
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        log_f = open(log_path, "w")
        proc  = subprocess.Popen(
            cmd, stdout=log_f, stderr=subprocess.STDOUT, text=True, env=env,
        )
    except FileNotFoundError:
        st.error(
            "Could not find the `opt` executable. "
            "Make sure you have run `pip install .` inside the `opt` conda environment "
            "and that the environment is active."
        )
        return

    threading.Thread(target=_wait_for_proc, args=(proc, rc_path), daemon=True).start()

    st.session_state["last_cmd"]              = " ".join(cmd)
    st.session_state["opt_running"]           = True
    st.session_state["opt_multi_mode"]        = False
    st.session_state["run_annotation_preset"] = module_inputs.get("annotation_preset", "")
    st.session_state["opt_log_path"]   = log_path
    st.session_state["opt_rc_path"]    = rc_path
    st.session_state["run_returncode"] = None
    st.session_state["run_stdout"]     = None
    st.session_state["out_dir"]        = out_dir
    st.session_state["run_module"]     = "all"
    st.session_state["current_step"]   = 0
    st.session_state["biotype_filter"] = "All"


def _run_multi_opt(global_args: dict, module_inputs: dict) -> None:
    base_out = global_args["out_dir"]
    os.makedirs(base_out, exist_ok=True)
    specs = build_multi_commands(global_args, module_inputs)
    for spec in specs:
        os.makedirs(spec["out_dir"], exist_ok=True)
        spec["log_path"] = os.path.join(spec["out_dir"], "_opt_run.log")
        spec["rc_path"]  = os.path.join(spec["out_dir"], "_opt_run.rc")
        for p in [spec["log_path"], spec["rc_path"]]:
            if os.path.exists(p):
                os.remove(p)

    first, *rest = specs
    try:
        _start_run(first)
    except FileNotFoundError:
        st.error(
            "Could not find the `opt` executable. "
            "Make sure you have run `pip install .` inside the `opt` conda environment."
        )
        return

    ref_names = [s["name"] for s in specs]
    st.session_state.update({
        "opt_running":        True,
        "opt_multi_mode":     True,
        "opt_current_run":    first,
        "opt_runs_queue":     rest,
        "opt_completed_runs": [],
        "opt_out_dir_base":   base_out,
        "out_dir":            base_out,
        "run_returncode":     None,
        "run_stdout":         None,
        "run_module":         "all",
        "current_step":       0,
        "biotype_filter":     "All",
        "last_cmd":           f"Multi-annotation: {ref_names}",
    })


def _render_running_status() -> None:
    """Dispatch to single or multi running status display."""
    if st.session_state.get("opt_multi_mode"):
        _render_multi_running_status()
    else:
        _render_single_running_status()


def _render_single_running_status() -> None:
    """Show live log while a single OPT run is in progress; auto-rerun until done."""
    log_path = st.session_state.get("opt_log_path", "")
    rc_path  = st.session_state.get("opt_rc_path", "")

    log_content = ""
    if log_path and os.path.exists(log_path):
        with open(log_path) as f:
            log_content = f.read()

    current_step = st.session_state.get("current_step", 0)
    if "FLIP module"  in log_content and current_step < 1:
        current_step = 1
    if "TRACK module" in log_content and current_step < 2:
        current_step = 2
    if "STAT module"  in log_content and current_step < 3:
        current_step = 3
    st.session_state["current_step"] = current_step

    if rc_path and os.path.exists(rc_path):
        with open(rc_path) as f:
            rc = int(f.read().strip())
        st.session_state["run_returncode"] = rc
        st.session_state["run_stdout"]     = log_content
        st.session_state["current_step"]   = 4 if rc == 0 else current_step
        st.session_state["opt_running"]    = False
        st.rerun()
        return

    with st.status("Running OPT...", expanded=True):
        st.markdown(_step_html(current_step), unsafe_allow_html=True)
        formatted = [_fmt_line(ln) for ln in log_content.splitlines()]
        formatted = [f for f in formatted if f]
        st.markdown("<br>".join(formatted), unsafe_allow_html=True)

    time.sleep(1)
    st.rerun()


def _render_multi_running_status() -> None:
    """Show multi-annotation progress; advance queue or merge when each run finishes."""
    completed = st.session_state.get("opt_completed_runs", [])
    current   = st.session_state.get("opt_current_run", {})
    queue     = st.session_state.get("opt_runs_queue", [])
    rc_path   = st.session_state.get("opt_rc_path", "")
    log_path  = st.session_state.get("opt_log_path", "")

    log_content = ""
    if log_path and os.path.exists(log_path):
        with open(log_path) as f:
            log_content = f.read()

    current_step = st.session_state.get("current_step", 0)
    for tag, n in [("FLIP module", 1), ("TRACK module", 2), ("STAT module", 3)]:
        if tag in log_content and current_step < n:
            current_step = n
    st.session_state["current_step"] = current_step

    if rc_path and os.path.exists(rc_path):
        with open(rc_path) as f:
            rc = int(f.read().strip())

        if rc != 0:
            st.session_state.update({
                "run_returncode": rc, "run_stdout": log_content,
                "opt_running": False, "opt_multi_mode": False,
            })
            st.rerun()
            return

        completed.append({**current, "log": log_content})
        st.session_state["opt_completed_runs"] = completed

        if queue:
            next_run = queue.pop(0)
            st.session_state["opt_runs_queue"]  = queue
            st.session_state["opt_current_run"] = next_run
            st.session_state["current_step"]    = 0
            _start_run(next_run)
            st.rerun()
            return
        else:
            base_out = st.session_state["opt_out_dir_base"]
            _merge_multi_results(completed, base_out)
            st.session_state.update({
                "run_returncode": 0,
                "run_stdout":     "Multi-annotation run complete.",
                "current_step":   4,
                "opt_running":    False,
                "opt_multi_mode": False,
                "biotype_filter": "All",
            })
            st.rerun()
            return

    # Still running — build progress header and display log
    done_names = {r["name"] for r in completed}
    parts = []
    for name in ["GENCODE", "CHESS", "RefSeq"]:
        color = _REF_COLORS[name]
        if name in done_names:
            parts.append(
                f'<span style="color:#2ca02c">&#10003; '
                f'<span style="background:{color};color:white;padding:1px 5px;'
                f'border-radius:3px;font-size:0.85em">{name}</span></span>')
        elif name == current.get("name"):
            parts.append(
                f'<span style="color:#1f77b4;font-weight:bold">&#9654; '
                f'<span style="background:{color};color:white;padding:1px 5px;'
                f'border-radius:3px;font-size:0.85em">{name}</span></span>')
        else:
            parts.append(
                f'<span style="color:#aaa">&#9675; '
                f'<span style="background:#aaa;color:white;padding:1px 5px;'
                f'border-radius:3px;font-size:0.85em">{name}</span></span>')
    progress_html = " &nbsp;&rarr;&nbsp; ".join(parts)

    with st.status(f"Running OPT — {current.get('name', '')} annotation...", expanded=True):
        st.markdown(progress_html, unsafe_allow_html=True)
        st.markdown(_step_html(current_step), unsafe_allow_html=True)
        formatted = [f for f in (_fmt_line(ln) for ln in log_content.splitlines()) if f]
        st.markdown("<br>".join(formatted), unsafe_allow_html=True)

    time.sleep(1)
    st.rerun()


# ── Results helpers ───────────────────────────────────────────────────────────
def _load_offtarget_data(out_dir: str):
    """Returns (summary_df, probes_df, all_probes_df). Any may be None if missing."""
    summary_path   = os.path.join(out_dir, "collapsed_summary_offtargets.tsv")
    probes_path    = os.path.join(out_dir, "probe2targets_offtargets.tsv")
    all_path       = os.path.join(out_dir, "probe2targets.tsv")
    summary_df     = pd.read_csv(summary_path, sep="\t") if os.path.exists(summary_path) else None
    probes_df      = pd.read_csv(probes_path,  sep="\t") if os.path.exists(probes_path)  else None
    all_probes_df  = pd.read_csv(all_path,     sep="\t") if os.path.exists(all_path)     else None
    return summary_df, probes_df, all_probes_df


def _count_file_lines(path: str) -> int:
    if not os.path.exists(path):
        return 0
    with open(path) as fh:
        return sum(1 for ln in fh if ln.strip())


def _build_ot_rows(summary_df, probes_df, all_probes_df=None) -> list:
    """Build one dict per (target_gene, off_target_gene) pair.

    Parses gene_names / transcript_types / cigars as parallel lists — columns
    present in both old (7-col) and new (11-col) probe2targets_offtargets.tsv.
    Off-target positions are those where gene_names[i] != probe_gene.

    Returns list of dicts with keys:
      target_gene, off_target_probes, off_target_gene, biotypes (set), cigars (set)
    """
    agg: dict = {}  # {target_gene: {ot_gene: {biotypes, cigars, probe_ids}}}

    if probes_df is not None:
        has_pg_col  = "probe_gene" in probes_df.columns
        has_ref_col = "reference_annotation" in probes_df.columns
        for _, row in probes_df.iterrows():
            if has_pg_col:
                probe_gene = str(row["probe_gene"])
            else:
                parts = str(row.get("probe_id", "")).split("|")
                probe_gene = parts[1] if len(parts) >= 2 else str(row.get("probe_id", ""))
            probe_id   = str(row.get("probe_id", ""))
            ref_name   = str(row.get("reference_annotation", "")).strip() if has_ref_col else ""
            gene_names = str(row.get("gene_names",       "")).strip("[]").split(",")
            ttypes     = str(row.get("transcript_types", "")).strip("[]").split(",")
            cigars_l   = str(row.get("cigars",           "")).strip("[]").split(",")

            for gname, ttype, cig in zip(gene_names, ttypes, cigars_l):
                gname = gname.strip()
                ttype = ttype.strip()
                cig   = cig.strip()
                if not gname or gname == probe_gene:
                    continue
                ot_entry = agg.setdefault(probe_gene, {}).setdefault(
                    gname, {"biotypes": set(), "cigars": set(), "probe_ids": set(),
                            "refs": set(), "ref_biotypes": {}}
                )
                ot_entry["biotypes"].add(ttype)
                ot_entry["cigars"].add(cig)
                ot_entry["probe_ids"].add(probe_id)
                if ref_name:
                    ot_entry["refs"].add(ref_name)
                    ot_entry["ref_biotypes"].setdefault(ref_name, set()).add(ttype)

    # Also scan probe2targets.tsv (all probes) for exclusive off-target hits:
    # probes that align ONLY to a non-target gene (n_genes=1, gene != probe_gene).
    # These don't appear in probe2targets_offtargets.tsv but are flagged in the summary.
    if all_probes_df is not None:
        has_ref_col_all = "reference_annotation" in all_probes_df.columns
        for _, row in all_probes_df.iterrows():
            parts = str(row.get("probe_id", "")).split("|")
            probe_gene = parts[1] if len(parts) >= 2 else str(row.get("probe_id", ""))
            probe_id   = str(row.get("probe_id", ""))
            ref_name   = str(row.get("reference_annotation", "")).strip() if has_ref_col_all else ""
            gene_names = str(row.get("gene_names",       "")).strip("[]").split(",")
            ttypes     = str(row.get("transcript_types", "")).strip("[]").split(",")
            cigars_l   = str(row.get("cigars",           "")).strip("[]").split(",")
            for gname, ttype, cig in zip(gene_names, ttypes, cigars_l):
                gname = gname.strip(); ttype = ttype.strip(); cig = cig.strip()
                if not gname or gname == probe_gene:
                    continue
                ot_entry = agg.setdefault(probe_gene, {}).setdefault(
                    gname, {"biotypes": set(), "cigars": set(), "probe_ids": set(),
                            "refs": set(), "ref_biotypes": {}}
                )
                ot_entry["biotypes"].add(ttype)
                ot_entry["cigars"].add(cig)
                ot_entry["probe_ids"].add(probe_id)
                if ref_name:
                    ot_entry["refs"].add(ref_name)
                    ot_entry["ref_biotypes"].setdefault(ref_name, set()).add(ttype)

    rows: list = []
    if summary_df is None:
        return rows

    for _, srow in summary_df.iterrows():
        tg      = str(srow["target_gene"])
        ot_data = agg.get(tg, {})

        if ot_data:
            for ot_gene, data in sorted(ot_data.items()):
                rows.append({
                    "target_gene":       tg,
                    "off_target_probes": len(data["probe_ids"]),
                    "off_target_gene":   ot_gene,
                    "biotypes":          data["biotypes"],
                    "cigars":            data["cigars"],
                    "refs":              data.get("refs", set()),
                    "ref_biotypes":      data.get("ref_biotypes", {}),
                })
        else:
            # Fallback: no probe-level data — derive off-target genes from summary row
            aligned_raw = str(srow.get("aligned_to", "")).strip("[]")
            for ot_gene in aligned_raw.split(","):
                ot_gene = ot_gene.strip()
                if ot_gene and ot_gene != tg:
                    rows.append({
                        "target_gene":       tg,
                        "off_target_probes": int(srow.get("n", 0)),
                        "off_target_gene":   ot_gene,
                        "biotypes":          set(),
                        "cigars":            set(),
                        "refs":              set(),
                        "ref_biotypes":      {},
                    })
    return rows


# ── Results dashboard ─────────────────────────────────────────────────────────
def render_results() -> None:
    if st.session_state["run_returncode"] is None:
        return

    st.divider()

    # Run metadata
    with st.expander("Run details", expanded=False):
        step_done = st.session_state.get("current_step", 0)
        st.markdown(_step_html(step_done), unsafe_allow_html=True)
        st.code(st.session_state.get("last_cmd", ""), language="bash")
        with st.expander("Full log", expanded=False):
            st.text(st.session_state.get("run_stdout") or "(no output)")

    out_dir = st.session_state["out_dir"]
    if not out_dir or not os.path.isdir(out_dir):
        st.warning("Output directory not found.")
        return

    if st.session_state["run_returncode"] != 0:
        st.error("OPT did not complete successfully. Check the log in 'Run details' above.")
        return

    summary_df, probes_df, all_probes_df = _load_offtarget_data(out_dir)
    if summary_df is None:
        st.info("No off-target summary file found in the output directory.")
        return

    # For single-annotation runs with a known preset, tag probes with the source name
    # so the Source and Gene biotype columns are populated correctly.
    # Read from disk so this works after page refresh or app restart.
    annotation_preset = st.session_state.get("run_annotation_preset", "")
    if not annotation_preset:
        preset_file = os.path.join(out_dir, "_opt_annotation_preset.txt")
        if os.path.exists(preset_file):
            with open(preset_file) as _f:
                annotation_preset = _f.read().strip()
            st.session_state["run_annotation_preset"] = annotation_preset
    if annotation_preset in _REF_COLORS:
        if probes_df is not None and "reference_annotation" not in probes_df.columns:
            probes_df = probes_df.copy()
            probes_df["reference_annotation"] = annotation_preset
        if all_probes_df is not None and "reference_annotation" not in all_probes_df.columns:
            all_probes_df = all_probes_df.copy()
            all_probes_df["reference_annotation"] = annotation_preset

    st.subheader("Predicted Off-Target Results")

    # Build the unified off-target row list (one entry per target→ot_gene pair)
    ot_rows = _build_ot_rows(summary_df, probes_df, all_probes_df)

    # ── Metric cards ──────────────────────────────────────────────────────────
    n_ot_genes  = len(summary_df)
    n_ot_probes = _count_file_lines(os.path.join(out_dir, "stat_off_target_probes.txt"))
    high_genes  = {r["target_gene"] for r in ot_rows
                   if any(bt in ("protein_coding", "mRNA") for bt in r["biotypes"])}
    n_high      = len(high_genes)

    mc1, mc2, mc3 = st.columns(3)
    mc1.metric("Genes with predicted off-target binding",                   n_ot_genes)
    mc2.metric("Genes with predicted protein-coding off-targets", n_high)
    mc3.metric("Probes with predicted off-target binding",                  n_ot_probes)

    st.divider()

    # ── Gene-level table ──────────────────────────────────────────────────────
    st.markdown("**Gene-level predicted off-target summary**")

    # Canonical groups for filter buttons (deduplicated, pseudogenes merged)
    canonical_groups: set = {_normalize_biotype(bt)
                             for r in ot_rows for bt in r["biotypes"] if bt}
    filter_options = ["All"] + sorted(canonical_groups)
    if st.session_state.get("biotype_filter", "All") not in filter_options:
        st.session_state["biotype_filter"] = "All"

    btn_cols = st.columns(len(filter_options))
    for i, grp in enumerate(filter_options):
        with btn_cols[i]:
            is_active = st.session_state.get("biotype_filter", "All") == grp
            btn_key = "filter_" + grp.replace(" ", "_").replace("/", "_")
            if st.button(grp, key=btn_key, use_container_width=True,
                         type="primary" if is_active else "secondary"):
                st.session_state["biotype_filter"] = grp
                st.rerun()

    active_filter = st.session_state.get("biotype_filter", "All")

    visible = [
        r for r in ot_rows
        if active_filter == "All"
        or any(_normalize_biotype(bt) == active_filter for bt in r["biotypes"])
    ]

    # ── Sort controls ─────────────────────────────────────────────────────────
    _SORT_COLS = ["Target gene", "Predicted off-target probes", "Predicted off-target genes", "Gene biotype"]
    _SORT_KEYS = {
        "Target gene":       lambda r: r["target_gene"].lower(),
        "Predicted off-target probes": lambda r: r["off_target_probes"],
        "Predicted off-target genes":   lambda r: r["off_target_gene"].lower(),
        "Gene biotype":      lambda r: min((_normalize_biotype(bt) for bt in r["biotypes"]), default=""),
    }
    sort_col_widget, sort_dir_widget = st.columns([3, 1])
    with sort_col_widget:
        sort_col = st.selectbox(
            "Sort by", _SORT_COLS,
            index=_SORT_COLS.index(st.session_state.get("table_sort_col", "Target gene")),
            key="table_sort_col",
        )
    with sort_dir_widget:
        st.write("")
        asc = st.toggle("Ascending", value=st.session_state.get("table_sort_asc", True),
                        key="table_sort_asc")

    visible = sorted(visible, key=_SORT_KEYS[sort_col], reverse=not asc)

    if not visible:
        st.info(f"No predicted off-target genes match the selected biotype filter: {active_filter.replace('_', ' ')}")
    else:
        th = ('background:#d8d8d8;color:#111111;padding:6px 8px;'
              'font-weight:bold;border-bottom:2px solid #bbb')
        html = (
            '<div style="background:#ffffff;border-radius:6px;overflow:auto;'
            'border:1px solid #ddd;padding:4px">'
            '<table style="width:100%;border-collapse:collapse;font-size:0.9em;color:#111111">'
            f'<thead><tr>'
            f'<th style="text-align:left;{th}">Target gene</th>'
            f'<th style="text-align:right;{th}">Predicted off-target probes</th>'
            f'<th style="text-align:left;{th}">Predicted off-target genes</th>'
            f'<th style="text-align:left;{th}">Gene biotype</th>'
            f'<th style="text-align:left;{th}">CIGAR</th>'
            f'<th style="text-align:left;{th}">Source</th>'
            f'</tr></thead><tbody>'
        )
        for i, r in enumerate(visible):
            bg = "#f5f5f5" if i % 2 == 0 else "#ffffff"
            # Build Gene biotype and Source columns
            ref_biotypes = r.get("ref_biotypes", {})
            if ref_biotypes:
                # Multi-annotation: pair each source with its biotype(s) in fixed order
                bt_parts, src_parts = [], []
                for ref in ["GENCODE", "CHESS", "RefSeq"]:
                    if ref not in ref_biotypes:
                        continue
                    color = _REF_COLORS.get(ref, "#999")
                    src_parts.append(
                        f'<span style="background:{color};color:white;padding:1px 5px;'
                        f'border-radius:3px;font-size:0.78em">{ref}</span>'
                    )
                    seen_n: set = set()
                    cell_badges = []
                    for bt in sorted(ref_biotypes[ref]):
                        norm = _normalize_biotype(bt)
                        if norm not in seen_n:
                            seen_n.add(norm)
                            cell_badges.append(_biotype_badge(bt))
                    # Stack multiple biotypes within a ref vertically; single is inline
                    bt_parts.append("<br>".join(cell_badges) if cell_badges else "—")
                _sep = ' <span style="color:#bbb;margin:0 4px">|</span> '
                bt_html  = _sep.join(bt_parts)
                ref_html = _sep.join(src_parts)
            else:
                # Single annotation: existing flat badge display
                seen_norms: set = set()
                bt_badges = []
                for bt in sorted(r["biotypes"]):
                    norm = _normalize_biotype(bt)
                    if norm not in seen_norms:
                        seen_norms.add(norm)
                        bt_badges.append(_biotype_badge(bt))
                bt_html = " ".join(bt_badges) if bt_badges else ""
                ref_html = " ".join(
                    f'<span style="background:{_REF_COLORS.get(ref, "#999")};color:white;'
                    f'padding:1px 5px;border-radius:3px;font-size:0.78em">{ref}</span>'
                    for ref in sorted(r.get("refs", set())) if ref
                )
            cigar_html = " ".join(
                f'<code style="background:#f0f0f0;color:#333;padding:1px 4px;'
                f'border-radius:3px;font-size:0.8em">{c}</code>'
                for c in sorted(r["cigars"])
            )
            html += (
                f'<tr style="background:{bg};border-bottom:1px solid #e8e8e8;color:#111111">'
                f'<td style="padding:5px 8px;font-weight:bold;color:#111111">{r["target_gene"]}</td>'
                f'<td style="padding:5px 8px;text-align:right;color:#111111">{r["off_target_probes"]}</td>'
                f'<td style="padding:5px 8px;color:#111111">{r["off_target_gene"]}</td>'
                f'<td style="padding:5px 8px">{bt_html}</td>'
                f'<td style="padding:5px 8px">{cigar_html}</td>'
                f'<td style="padding:5px 8px">{ref_html}</td>'
                f'</tr>'
            )
        html += '</tbody></table></div>'
        st.markdown(html, unsafe_allow_html=True)
        st.caption(f"Showing {len(visible)} rows ({n_ot_genes} target genes with off-target binding).")

    summary_path = os.path.join(out_dir, "collapsed_summary_offtargets.tsv")
    if os.path.exists(summary_path):
        with open(summary_path, "rb") as fh:
            st.download_button(
                "Download gene-level summary (collapsed_summary_offtargets.tsv)",
                data=fh, file_name="collapsed_summary_offtargets.tsv",
                key="dl_collapsed_ot",
            )

    st.divider()

    # ── Probe-level detail ────────────────────────────────────────────────────
    with st.expander("Probe-level predicted off-target detail", expanded=False):
        probes_path = os.path.join(out_dir, "probe2targets_offtargets.tsv")
        if probes_df is not None:
            # Rebuild off-target gene names, types, and CIGARs from the same
            # parallel lists so all three columns always have the same count.
            # Deduplicates by (gene, type, cigar) triple; uses | as delimiter.
            has_pg = "probe_gene" in probes_df.columns

            def _ot_info(row):
                probe_gene = str(row["probe_gene"]) if has_pg else (
                    str(row.get("probe_id", "")).split("|")[1]
                    if len(str(row.get("probe_id", "")).split("|")) >= 2
                    else str(row.get("probe_id", ""))
                )
                genes  = str(row.get("gene_names",       "")).strip("[]").split(",")
                ttypes = str(row.get("transcript_types",  "")).strip("[]").split(",")
                cigars = str(row.get("cigars",            "")).strip("[]").split(",")
                seen = set()
                ot_genes, ot_types, ot_cigars = [], [], []
                for g, t, c in zip(genes, ttypes, cigars):
                    g, t, c = g.strip(), t.strip(), c.strip()
                    if not g or g == probe_gene:
                        continue
                    key = (g, t, c)
                    if key not in seen:
                        seen.add(key)
                        ot_genes.append(g)
                        ot_types.append(t)
                        ot_cigars.append(c)
                return pd.Series({
                    "_ot_genes":  " | ".join(ot_genes),
                    "_ot_types":  " | ".join(ot_types),
                    "_ot_cigars": " | ".join(ot_cigars),
                })

            display_df = probes_df.copy()
            ot_cols = display_df.apply(_ot_info, axis=1)
            display_df["_ot_genes"]  = ot_cols["_ot_genes"]
            display_df["_ot_types"]  = ot_cols["_ot_types"]
            display_df["_ot_cigars"] = ot_cols["_ot_cigars"]

            # Select and rename columns for display
            col_map = {
                "probe_id":             "Probe ID",
                "probe_gene":           "Probe gene",
                "_ot_genes":            "Off-target genes",
                "_ot_types":            "Off-target biotypes",
                "_ot_cigars":           "Off-target CIGARs",
                "reference_annotation": "Source",
            }
            show_cols = [c for c in col_map if c in display_df.columns]
            display_df = display_df[show_cols].rename(columns=col_map)

            total = len(display_df)
            max_rows = st.slider(
                "Rows to show", min_value=10, max_value=max(total, 10),
                value=min(50, total), step=10, key="probe_detail_slider",
            )
            st.dataframe(display_df.head(max_rows), use_container_width=True)
            st.caption(f"Showing {min(max_rows, total)} of {total} probes.")
        else:
            st.info("probe2targets_offtargets.tsv not found in output directory.")
        if os.path.exists(probes_path):
            with open(probes_path, "rb") as fh:
                st.download_button(
                    "Download probe-level predicted off-targets (probe2targets_offtargets.tsv)",
                    data=fh, file_name="probe2targets_offtargets.tsv",
                    key="dl_probes_ot",
                )

    st.divider()

    # ── Output file descriptions ───────────────────────────────────────────────
    st.subheader("Output File Descriptions")
    st.caption(f"All output files are saved to: `{out_dir}`")
    st.markdown(
        """
| File | Description |
|---|---|
| `collapsed_summary_offtargets.tsv` | One row per target gene with off-target binding — lists the off-target genes, probe counts, and biotypes. |
| `probe2targets_offtargets.tsv` | One row per probe that maps to more than one gene — includes all alignment targets, CIGAR strings, and biotypes. |
| `collapsed_summary.tsv` | Summary across all target genes (including those with no off-target binding). |
| `probe2targets.tsv` | Full alignment table for all probes, including on-target hits. |
| `stat_off_target_probes.txt` | List of probe IDs with at least one off-target alignment. |
| `stat_off_target_genes.txt` | List of off-target gene names detected across all probes. |
| `stat_missed_probes.txt` | Probes that did not align to their intended target gene. |
| `stat_missed_genes.txt` | Target genes for which no probes aligned. |
| `fwd_oriented.fa` | Strand-corrected probe sequences produced by the flip step. |
        """,
        unsafe_allow_html=False,
    )


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    init_session_state()
    render_header()
    st.divider()

    global_args, module_inputs = render_all_inputs()
    st.divider()

    run_col, load_col = st.columns([3, 1])
    with run_col:
        if st.button("Run OPT", type="primary", use_container_width=True):
            run_opt(global_args, module_inputs)
    with load_col:
        if st.button("Load previous results", use_container_width=True,
                     help="Load results from a previous run without re-running OPT. "
                          "Set the Output Directory above to the folder from your previous run, "
                          "then click this button."):
            out_dir = st.session_state.get("all_out_dir_val", "./opt_results").strip()
            summary = os.path.join(out_dir, "collapsed_summary_offtargets.tsv")
            if os.path.isdir(out_dir) and os.path.exists(summary):
                st.session_state["run_returncode"] = 0
                st.session_state["out_dir"]        = out_dir
                st.session_state["run_stdout"]     = ""
                st.session_state["last_cmd"]       = "(loaded from previous run)"
                st.session_state["current_step"]   = 4
                st.session_state["biotype_filter"] = "All"
            else:
                st.warning("No results found in the output directory.")

    if st.session_state.get("opt_running"):
        _render_running_status()
    else:
        render_results()


if __name__ == "__main__":
    main()
