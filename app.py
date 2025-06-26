# main.py  ──────────────────────────────────────────────────────────────
# Extended version of the original app.py
# • Adds two top-level sections:  H Data  and  F Data
# • Each section contains Spleen, Kidney, Liver sub-tabs
# • Uses lazy construction to avoid loading data until a tab is opened

import os
import param
import requests
import panel as pn
from utils import load_and_prepare_data
from volcano_plot import generate_plot

pn.config.raw_css.append("""
html, body {
    margin: 0;
    padding: 0;
    overflow: hidden !important;
    height: 100%;
}

/* Applies to all tab buttons */
.bk-tab {
    font-size: 12px !important;
    font-weight: bold !important;
    padding: 8px 12px !important;
}
""")

# ────────────────────────────── Panel setup ───────────────────────────
pn.extension("plotly", "tabulator")  # single extension call keeps memory tight

ENV_CHECK = os.environ.get("ENV_CHECK", "")          # "DEV" for local dev
IFTTT_KEY = os.environ.get("IFTTT_API_KEY", "")      # optional deploy ping

# ────────────────────────── Core interactive class ────────────────────
class VolcanoApp(param.Parameterized):
    """Interactive compound explorer + volcano plot."""
    comparison = param.Integer(0)

    def __init__(self, data_file: str, distance_file: str, mapping_key: str, **params):
        super().__init__(**params)

        # Lazy-load and preprocess data for the chosen organ/dataset group
        self.df, self.comparisons = load_and_prepare_data(
            data_file, distance_file, mapping_key
        )

        # Build widgets
        self.comparison_names = [
            comp.get("title", f"Comparison {i+1}") for i, comp in enumerate(self.comparisons)
        ]
        self.comparison_select = pn.widgets.Select(
            name="Comparison",
            options=dict(zip(self.comparison_names, range(len(self.comparison_names)))),
        )
        self.search_input = pn.widgets.TextInput(placeholder="Search m/z…", name="Search")

        self.clear_all_button = pn.widgets.Button(name="Checkbox Clear", button_type="danger")
        self.reset_button = pn.widgets.Button(
            name="Reload Page",
            button_type="danger",
            width=100,
            height=30,
            styles={
                "background": "red",
                "color": "white",
                "font-size": "8px",
                "font-weight": "bold",
                "border-radius": "4px",
            },
        )

        self.table_columns = [
            {"field": "m/z", "title": "m/z", "headerTooltip": "m/z value"},
            {"field": "RT [min]", "title": "RT [min]", "headerTooltip": "Retention time"},
            {"field": "Formula", "title": "Formula", "headerTooltip": "Chemical formula"},
        ]
        self.table = pn.widgets.Tabulator(
            pagination=None,
            #height=725,
            selectable="checkbox",
            header_align="center",
            layout="fit_columns",
            sizing_mode="stretch_both",
            show_index=False,
            theme="midnight",
            hidden_columns=["Compounds ID", "abs_fc"],
        )
        self.table.columns = self.table_columns

        self.plot_pane = pn.pane.Plotly(sizing_mode="stretch_both") #, min_height=900)

        # Wire callbacks
        self.comparison_select.param.watch(self._update_comparison, "value")
        self.table.param.watch(self._update_selection, "selection")
        self.search_input.param.watch(self._apply_filter, "value")
        self.clear_all_button.on_click(self._clear_all)
        self.reset_button.on_click(self._reset_app)

        # First render
        self._update_comparison(None)

    # ─────────── Callbacks ───────────
    def _update_comparison(self, _):
        idx = self.comparison_select.value
        self.comparison = idx

        fc_col = self.comparisons[idx]["fold_change_col"]
        pv_col = self.comparisons[idx]["p_value_col"]

        table_data = (
            self.df[["Compounds ID", "m/z", "RT [min]", "Formula", fc_col, pv_col]]
            .dropna(subset=[fc_col])
            .copy()
            .rename(
                columns={
                    "RT [min]": "RT",
                    pv_col: "P-val",
                    fc_col: "Log2 Fold",
                }
            )
        )
        table_data["abs_fc"] = table_data["Log2 Fold"].abs()
        table_data = table_data.sort_values("abs_fc", ascending=False)

        self.df_filtered = table_data            # cache for search
        self.search_input.value = ""             # clear search bar
        self.table.value = table_data
        self.table.selection = []                # clear any checks

        self.plot_pane.object = generate_plot(self.df, idx, self.comparisons)

        # Update visible columns once (keeps header tool-tips)
        self.table.columns = [
            {"field": "m/z", "title": "m/z", "headerTooltip": "m/z value"},
            {"field": "RT", "title": "RT [min]", "headerTooltip": "Retention time"},
            {"field": "Formula", "title": "Formula", "headerTooltip": "Chemical formula"},
            {"field": "P-val", "title": "P-val", "headerTooltip": "P-value"},
            {"field": "Log2 Fold", "title": "Log2 Fold", "headerTooltip": "Log2 fold-change"},
        ]

    def _apply_filter(self, event):
        query = event.new.strip()
        self.table.value = (
            self.df_filtered if not query else self.df_filtered[self.df_filtered["m/z"].astype(str).str.contains(query)]
        )

    def _update_selection(self, _):
        selected_ids = []
        if self.table.selection:
            sel_rows = [self.table.value.iloc[i] for i in self.table.selection]
            selected_ids = [str(row["Compounds ID"]) for row in sel_rows]

        self.plot_pane.object = generate_plot(
            self.df, self.comparison, self.comparisons, selected_ids
        )

    def _clear_all(self, _):        # quick de-select button
        self.table.selection = []

    def _reset_app(self, _):        # reload current organ
        self.search_input.value = ""
        self.comparison_select.value = 0
        self.table.selection = []
        self._update_comparison(None)

    # ─────────── Layout builder ───────────
    def panel(self):
        button_row = pn.Row(
            self.clear_all_button, self.reset_button, sizing_mode="stretch_width"
        )

        control_panel = pn.Column(
            "## Compound Explorer",
            self.comparison_select,
            self.search_input,
            button_row,
            self.table,
            styles={"background": "#606060"},
            width=600,
            margin=0,
        )

        main_panel = pn.Column(
            "## Volcano Plot",
            self.plot_pane,
            styles={"background": "#606060"},
            min_width=550,
            margin=0,
            sizing_mode="stretch_width",
        )

        return pn.Row(
            control_panel,
            main_panel,
            styles={"background": "#606060", "height": "100vh"},
            sizing_mode="stretch_width",
            margin=0,
        )

# ───────────────────────────── Deploy helper ───────────────────────────
def notify_webhook() -> None:
    """Optional IFTTT ping when deploying in prod."""
    if ENV_CHECK != "DEV" and IFTTT_KEY:
        webhook_url = f"https://maker.ifttt.com/trigger/sidebar/json/with/key/{IFTTT_KEY}"
        try:
            requests.post(webhook_url, json={"value1": "Deploying Website"}, timeout=4).raise_for_status()
        except requests.exceptions.RequestException as exc:
            # Print but don’t block startup
            print(f"[WARN] IFTTT webhook failed: {exc}")

# ────────────────────────── Helper for tab groups ──────────────────────
def _organ_tabs(prefix: str) -> pn.Tabs:
    """
    Build inner Tabs (Spleen/Kidney/Liver) for either H Data or F Data.

    Parameters
    ----------
    prefix : str
        ""Re"  → ReSpleen.csv, ReSpleen_by_distance_named.csv, …
        "F_Re" → F_ReSpleen.csv, F_ReSpleen_by_distance_named.csv, …
    """
    make_panel = lambda organ: pn.panel(
        lambda: VolcanoApp(
            f"{prefix}{organ}.csv",
            f"{prefix}{organ}_by_distance_named.csv",
            f"{prefix}{organ}.csv",
        ).panel()
    )

    return pn.Tabs(
        ("Spleen", make_panel("Spleen")),
        ("Kidney", make_panel("Kidney")),
        ("Liver",  make_panel("Liver")),
        dynamic=True,
    )

def _single_tab(prefix: str, label: str) -> pn.Tabs:
    """
    Builds a single tab (no organs) for datasets like Plasma.
    """
    panel = pn.panel(
        lambda: VolcanoApp(
            f"{prefix}.csv",
            f"{prefix}_by_distance_named.csv",
            f"{prefix}.csv",
        ).panel()
    )
    return pn.Tabs((label, panel), dynamic=True)

# ──────────────────────────────── main() ───────────────────────────────
def main() -> None:
    """Serve the Panel application on the desired port."""
    root_tabs = pn.Tabs(
        ("H Data", _organ_tabs("Re")),
        ("F Data", _organ_tabs("F_Re")),
        ("Plasma Data", _single_tab("SysPlas2", "Plasma")),
        dynamic=True,
        sizing_mode="stretch_height",
    )

    port = 4603 if ENV_CHECK == "DEV" else 80
    pn.serve(root_tabs, port=port, websocket_origin=["*"])

# ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    notify_webhook()
    main()