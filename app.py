import streamlit as st
import subprocess
import os
import re
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


# ── Biotype color helpers ─────────────────────────────────────────────────────
def _biotype_color(biotype: str) -> str:
    if biotype in ("protein_coding", "mRNA"):
        return "#d62728"
    if "pseudogene" in biotype:
        return "#888888"
    if biotype == "lncRNA":
        return "#e07b00"
    if biotype in ("miscRNA", "miRNA", "snoRNA", "snRNA"):
        return "#9467bd"
    return "#1f77b4"


def _biotype_badge(biotype: str) -> str:
    color = _biotype_color(biotype)
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
        # Path fields managed via Browse (not owned by any widget via key=)
        "all_query_path_val": "",
        "all_target_val":     "",
        "all_annotation_val": "",
        "all_syn_val":        "",
        "all_out_dir_val":    "./opt_results",
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


# ── Header ────────────────────────────────────────────────────────────────────
def render_header() -> None:
    st.title("OPT — Off-target Probe Tracker")
    st.markdown(
        """
        OPT identifies off-target binding of spatial transcriptomics probes
        (Xenium, MERSCOPE, CosMx, and others) against a reference transcriptome
        using nucleotide alignment (nucmer). The pipeline runs three steps automatically:
        **flip** corrects probe strand orientation by aligning to the source transcriptome,
        **track** aligns probes to all transcripts in the target transcriptome and flags
        probes that map to more than one gene as off-target, and **stat** summarizes
        off-target binding by gene and transcript biotype.

        **Citation:** Hallinan et al., *eLife* 2025.
        [https://elifesciences.org/reviewed-preprints/107070](https://elifesciences.org/reviewed-preprints/107070)
        """
    )
    st.info(
        "**Probe FASTA header format required:** `>gene_id|gene_name|accession`  \n"
        "Example: `>ENSG00000170458|CD14|22f9405`"
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
        st.write("")  # vertical alignment spacer
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
        threads = st.slider("Threads", min_value=1, max_value=16, value=1)

    st.divider()
    st.subheader("Input Files")

    query_path = _path_field(
        "Probe FASTA", "all_query_path_val", "all_query_browse",
        "Select Probe FASTA",
        "/path/to/probes.fa or probes.fa.gz",
        "Path to probe sequences (.fa or .fa.gz). Header format: >gene_id|gene_name|accession",
    )
    target = _path_field(
        "Target transcript FASTA", "all_target_val", "all_target_browse",
        "Select Target Transcript FASTA",
        "/path/to/transcripts.fa.gz",
        "Local path to target transcript sequences. Supports .fa and .fa.gz.",
    )
    annotation = _path_field(
        "Annotation GFF/GTF", "all_annotation_val", "all_annotation_browse",
        "Select Annotation GFF/GTF",
        "/path/to/annotation.gff.gz",
        "Local path to transcript annotation. Supports .gff, .gtf, and .gz.",
    )
    syn_file = _path_field(
        "Gene synonyms CSV", "all_syn_val", "all_syn_browse",
        "Select Gene Synonyms CSV",
        "/path/to/synonyms.csv",
        "Optional. Two-column CSV mapping gene name variants (e.g. MYCN,N-Myc).",
        optional=True,
    )

    st.divider()
    st.subheader("Analysis Options")

    col_a, col_b = st.columns(2)
    with col_a:
        pad_length = st.number_input(
            "Pad length (-pl)", min_value=0, value=0, step=1,
            help="Length of probe ends where misalignment is tolerated.",
        )
    with col_b:
        excl_pseudo = st.checkbox("Exclude pseudogenes (--exclude-pseudo)", value=False)
        pc_only     = st.checkbox("Protein-coding only (--pc-only)", value=False)

    with st.expander("Advanced Options", expanded=False):
        bam = st.checkbox(
            "Store alignments as BAM (--bam)", value=False,
            help="Save nucmer alignments as BAM instead of SAM.",
        )
        min_exact = st.number_input(
            "Min exact match length (-l)", min_value=1, value=20, step=1,
            help="Minimum exact match length for nucmer.",
        )
        keep_dot = st.checkbox(
            "Keep version numbers in gene IDs (--keep-dot)", value=False,
        )
        force = st.checkbox(
            "Force rebuild (--force)", value=False,
            help="Ignore cached results from previous runs and recompute everything.",
        )

        st.caption("GFF/GTF Schema (--schema)")
        st.caption(
            "Five fields used to parse the annotation file. "
            "Defaults work for GENCODE GFF. Change if using RefSeq or CHESS."
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

        st.info(
            "For Bowtie2 or a custom aligner, use the CLI version:  \n"
            "https://github.com/JEFworks-Lab/off-target-probe-tracker"
        )

    global_args = {
        "out_dir":       out_dir,
        "threads":       threads,
        "bam":           bam,
        "min_exact":     min_exact,
        "keep_dot":      keep_dot,
        "force":         force,
        "schema":        ",".join(schema_fields),
        "schema_fields": schema_fields,
    }
    module_inputs = {
        "query":          query_path,
        "target":         target,
        "annotation":     annotation,
        "syn_file":       syn_file,
        "pad_length":     pad_length,
        "exclude_pseudo": excl_pseudo,
        "pc_only":        pc_only,
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
    if global_args["keep_dot"]:
        cmd.append("--keep-dot")
    if global_args["force"]:
        cmd.append("--force")

    cmd.append("all")

    cmd += ["-q", module_inputs["query"].strip()]
    cmd += ["-t", module_inputs["target"].strip()]
    cmd += ["-a", module_inputs["annotation"].strip()]
    cmd += ["-pl", str(module_inputs["pad_length"])]

    if module_inputs.get("exclude_pseudo"):
        cmd.append("--exclude-pseudo")
    if module_inputs.get("pc_only"):
        cmd.append("--pc-only")
    syn = module_inputs.get("syn_file", "").strip()
    if syn:
        cmd += ["-s", syn]

    return cmd


# ── Validate inputs ───────────────────────────────────────────────────────────
def validate_inputs(global_args: dict, module_inputs: dict) -> list:
    errors = []
    if not global_args["out_dir"]:
        errors.append("Output directory is required.")
    for field in global_args.get("schema_fields", []):
        if "," in field:
            errors.append(f"Schema field contains a comma: '{field}'. Remove it.")
    for field, label in [
        ("query",      "Probe FASTA"),
        ("target",     "Target transcript FASTA"),
        ("annotation", "Annotation GFF/GTF"),
    ]:
        path = module_inputs.get(field, "").strip()
        if not path:
            errors.append(f"{label} path is required.")
        elif not os.path.exists(path):
            errors.append(f"{label} path does not exist: '{path}'")
    syn = module_inputs.get("syn_file", "").strip()
    if syn and not os.path.exists(syn):
        errors.append(f"Gene synonyms file does not exist: '{syn}'")
    if module_inputs.get("pc_only") and module_inputs.get("exclude_pseudo"):
        errors.append("Cannot use both --pc-only and --exclude-pseudo at the same time.")
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
def run_opt(global_args: dict, module_inputs: dict) -> None:
    errors = validate_inputs(global_args, module_inputs)
    if errors:
        for e in errors:
            st.error(e)
        return

    cmd = build_command(global_args, module_inputs)
    st.session_state["last_cmd"]     = " ".join(cmd)
    st.session_state["current_step"] = 0
    os.makedirs(global_args["out_dir"], exist_ok=True)

    collected_lines = []
    try:
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=env,
        )
    except FileNotFoundError:
        st.error(
            "Could not find the `opt` executable. "
            "Make sure you have run `pip install .` inside the `opt` conda environment "
            "and that the environment is active."
        )
        return

    with st.status("Running OPT...", expanded=True) as status:
        step_placeholder = st.empty()
        log_placeholder  = st.empty()
        formatted_lines  = []
        current_step     = 0

        for raw_line in proc.stdout:
            collected_lines.append(raw_line)
            clean = _strip_ansi(raw_line)
            # Detect which module just started
            if "FLIP module" in clean and "START" in clean:
                current_step = 1
            elif "TRACK module" in clean and "START" in clean:
                current_step = 2
            elif "STAT module" in clean and "START" in clean:
                current_step = 3
            step_placeholder.markdown(_step_html(current_step), unsafe_allow_html=True)
            fmt = _fmt_line(raw_line)
            if fmt:
                formatted_lines.append(fmt)
                log_placeholder.markdown("<br>".join(formatted_lines), unsafe_allow_html=True)

        proc.wait()
        rc = proc.returncode
        if rc == 0:
            current_step = 4
            step_placeholder.markdown(_step_html(4), unsafe_allow_html=True)
            status.update(label="OPT finished successfully", state="complete", expanded=False)
        else:
            status.update(label=f"OPT failed (exit code {rc})", state="error", expanded=True)

    st.session_state["run_stdout"]     = "".join(collected_lines)
    st.session_state["run_returncode"] = rc
    st.session_state["run_module"]     = "all"
    st.session_state["out_dir"]        = global_args["out_dir"]
    st.session_state["current_step"]   = current_step
    st.session_state["biotype_filter"] = "All"


# ── Results helpers ───────────────────────────────────────────────────────────
def _load_offtarget_data(out_dir: str):
    """Returns (summary_df, probes_df). Either may be None if file is missing."""
    summary_path = os.path.join(out_dir, "collapsed_summary_offtargets.tsv")
    probes_path  = os.path.join(out_dir, "probe2targets_offtargets.tsv")
    summary_df = pd.read_csv(summary_path, sep="\t") if os.path.exists(summary_path) else None
    probes_df  = pd.read_csv(probes_path,  sep="\t") if os.path.exists(probes_path)  else None
    return summary_df, probes_df


def _count_file_lines(path: str) -> int:
    if not os.path.exists(path):
        return 0
    with open(path) as fh:
        return sum(1 for ln in fh if ln.strip())


def _build_ot_rows(summary_df, probes_df) -> list:
    """Build one dict per (target_gene, off_target_gene) pair.

    Parses gene_names / transcript_types / cigars as parallel lists — columns
    present in both old (7-col) and new (11-col) probe2targets_offtargets.tsv.
    Off-target positions are those where gene_names[i] != probe_gene.

    Returns list of dicts with keys:
      target_gene, off_target_probes, off_target_gene, biotypes (set), cigars (set)
    """
    agg: dict = {}  # {target_gene: {ot_gene: {biotypes, cigars, probe_ids}}}

    if probes_df is not None:
        has_pg_col = "probe_gene" in probes_df.columns
        for _, row in probes_df.iterrows():
            if has_pg_col:
                probe_gene = str(row["probe_gene"])
            else:
                parts = str(row.get("probe_id", "")).split("|")
                probe_gene = parts[1] if len(parts) >= 2 else str(row.get("probe_id", ""))
            probe_id   = str(row.get("probe_id", ""))
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
                    gname, {"biotypes": set(), "cigars": set(), "probe_ids": set()}
                )
                ot_entry["biotypes"].add(ttype)
                ot_entry["cigars"].add(cig)
                ot_entry["probe_ids"].add(probe_id)

    rows: list = []
    if summary_df is None:
        return rows

    for _, srow in summary_df.iterrows():
        tg      = str(srow["target_gene"])
        ot_data = agg.get(tg, {})

        if ot_data:
            for ot_gene, data in sorted(ot_data.items()):
                rows.append({
                    "target_gene":      tg,
                    "off_target_probes": len(data["probe_ids"]),
                    "off_target_gene":  ot_gene,
                    "biotypes":         data["biotypes"],
                    "cigars":           data["cigars"],
                })
        else:
            # Fallback: no probe-level data — derive off-target genes from summary row
            aligned_raw = str(srow.get("aligned_to", "")).strip("[]")
            for ot_gene in aligned_raw.split(","):
                ot_gene = ot_gene.strip()
                if ot_gene and ot_gene != tg:
                    rows.append({
                        "target_gene":      tg,
                        "off_target_probes": int(srow.get("n", 0)),
                        "off_target_gene":  ot_gene,
                        "biotypes":         set(),
                        "cigars":           set(),
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

    summary_df, probes_df = _load_offtarget_data(out_dir)
    if summary_df is None:
        st.info("No off-target summary file found in the output directory.")
        return

    st.subheader("Off-Target Results")

    # Build the unified off-target row list (one entry per target→ot_gene pair)
    ot_rows = _build_ot_rows(summary_df, probes_df)

    # ── Metric cards ──────────────────────────────────────────────────────────
    n_ot_genes  = len(summary_df)
    n_ot_probes = _count_file_lines(os.path.join(out_dir, "stat_off_target_probes.txt"))
    high_genes  = {r["target_gene"] for r in ot_rows
                   if any(bt in ("protein_coding", "mRNA") for bt in r["biotypes"])}
    n_high      = len(high_genes)

    mc1, mc2, mc3 = st.columns(3)
    mc1.metric("Genes with off-target binding",                   n_ot_genes)
    mc2.metric("Probes with off-target binding",                  n_ot_probes)
    mc3.metric("Protein-coding off-targets", n_high)

    st.divider()

    # ── Gene-level table ──────────────────────────────────────────────────────
    st.markdown("**Gene-level off-target summary**")

    all_biotypes: set = {bt for r in ot_rows for bt in r["biotypes"]}
    all_biotypes.discard("")

    filter_options = ["All"] + sorted(all_biotypes)
    if st.session_state.get("biotype_filter", "All") not in filter_options:
        st.session_state["biotype_filter"] = "All"

    btn_cols = st.columns(len(filter_options))
    for i, bt in enumerate(filter_options):
        with btn_cols[i]:
            is_active = st.session_state.get("biotype_filter", "All") == bt
            label = bt if bt == "All" else bt.replace("_", " ")
            if st.button(label, key=f"filter_{bt}", use_container_width=True,
                         type="primary" if is_active else "secondary"):
                st.session_state["biotype_filter"] = bt

    active_filter = st.session_state.get("biotype_filter", "All")

    visible = [
        r for r in ot_rows
        if active_filter == "All" or active_filter in r["biotypes"]
    ]

    # ── Sort controls ─────────────────────────────────────────────────────────
    _SORT_COLS = ["Target gene", "Off-target probes", "Off-target gene", "Gene biotype"]
    _SORT_KEYS = {
        "Target gene":       lambda r: r["target_gene"].lower(),
        "Off-target probes": lambda r: r["off_target_probes"],
        "Off-target gene":   lambda r: r["off_target_gene"].lower(),
        "Gene biotype":      lambda r: min(r["biotypes"]) if r["biotypes"] else "",
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
        st.info(f"No off-target genes match the selected biotype filter: {active_filter.replace('_', ' ')}")
    else:
        th = ('background:#d8d8d8;color:#111111;padding:6px 8px;'
              'font-weight:bold;border-bottom:2px solid #bbb')
        html = (
            '<div style="background:#ffffff;border-radius:6px;overflow:auto;'
            'border:1px solid #ddd;padding:4px">'
            '<table style="width:100%;border-collapse:collapse;font-size:0.9em;color:#111111">'
            f'<thead><tr>'
            f'<th style="text-align:left;{th}">Target gene</th>'
            f'<th style="text-align:right;{th}">Off-target probes</th>'
            f'<th style="text-align:left;{th}">Off-target gene</th>'
            f'<th style="text-align:left;{th}">Gene biotype</th>'
            f'<th style="text-align:left;{th}">CIGAR</th>'
            f'</tr></thead><tbody>'
        )
        for i, r in enumerate(visible):
            bg = "#f5f5f5" if i % 2 == 0 else "#ffffff"
            bt_html = (
                " ".join(_biotype_badge(bt) for bt in sorted(r["biotypes"]))
                if r["biotypes"] else ""
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
    with st.expander("Probe-level off-target detail", expanded=False):
        probes_path = os.path.join(out_dir, "probe2targets_offtargets.tsv")
        if probes_df is not None:
            show_cols = [c for c in
                         ["probe_id", "probe_gene", "offtarget_gene_names",
                          "offtarget_gene_types", "concern_level"]
                         if c in probes_df.columns]
            if not show_cols:
                show_cols = list(probes_df.columns)
            total = len(probes_df)
            max_rows = st.slider(
                "Rows to show", min_value=10, max_value=max(total, 10),
                value=min(50, total), step=10, key="probe_detail_slider",
            )
            st.dataframe(probes_df[show_cols].head(max_rows), use_container_width=True)
            st.caption(f"Showing {min(max_rows, total)} of {total} probes.")
        else:
            st.info("probe2targets_offtargets.tsv not found in output directory.")
        if os.path.exists(probes_path):
            with open(probes_path, "rb") as fh:
                st.download_button(
                    "Download probe-level off-targets (probe2targets_offtargets.tsv)",
                    data=fh, file_name="probe2targets_offtargets.tsv",
                    key="dl_probes_ot",
                )

    st.divider()

    # ── Download all key output files ─────────────────────────────────────────
    st.markdown("**Download output files**")
    dl_files = {
        "collapsed_summary.tsv":     "All-gene summary",
        "probe2targets.tsv":         "All-probe alignments",
        "stat_missed_probes.txt":    "Missed-target probes",
        "stat_off_target_genes.txt": "Off-target gene list",
        "stat_missed_genes.txt":     "Missed-target gene list",
    }
    dl_cols = st.columns(len(dl_files))
    for col, (fname, label) in zip(dl_cols, dl_files.items()):
        fpath = os.path.join(out_dir, fname)
        with col:
            if os.path.exists(fpath):
                with open(fpath, "rb") as fh:
                    st.download_button(label, data=fh, file_name=fname, key=f"dl_{fname}")
            else:
                st.caption(f"*{fname}* not found")


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
        if st.button("Load previous results", use_container_width=True):
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

    render_results()


if __name__ == "__main__":
    main()
