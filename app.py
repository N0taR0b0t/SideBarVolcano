# app.py - Simplified layout with better proportions
import os
import param
import requests
import panel as pn
from utils import load_and_prepare_data
from volcano_plot import generate_plot

# Initialize Panel extension
pn.extension('plotly', 'tabulator')
ENV_CHECK = os.environ.get("ENV_CHECK", "")
IFTTT_KEY = os.environ.get("IFTTT_API_KEY", "")

class VolcanoApp(param.Parameterized):
    comparison = param.Integer(0)

    def __init__(self, data_file, distance_file, mapping_key, **params):
        super().__init__(**params)
        self.df, self.comparisons = load_and_prepare_data(data_file, distance_file, mapping_key)


        # Load data
        #self.df, self.comparisons = load_and_prepare_data()

        # Create comparison dropdown options
        self.comparison_names = [comp.get("title", f"Comparison {i+1}")
                                 for i, comp in enumerate(self.comparisons)]

        # Create widgets
        self.comparison_select = pn.widgets.Select(
            name='Comparison',
            options=dict(zip(self.comparison_names, range(len(self.comparison_names))))
        )

        # Create search bar
        self.search_input = pn.widgets.TextInput(
            placeholder='Search m/z...', name='Search'
        )

        # Create buttons
        #self.select_all_button = pn.widgets.Button(name="Select All", button_type="primary")
        self.clear_all_button = pn.widgets.Button(name="Checkbox Clear", button_type="danger")
        self.reset_button = pn.widgets.Button(
            name="Reload Page", 
            button_type="danger",
            width=100,
            height=30,
            styles={
                'background': 'red',
                'color': 'white',
                'font-size': '8px',
                'font-weight': 'bold',
                'border-radius': '4px'
            }
        )

        # Define columns for the table
        self.table_columns = [
            {'field': 'm/z', 'title': 'm/z', 'headerTooltip': 'm/z value'},
            {'field': 'RT [min]', 'title': 'RT [min]', 'headerTooltip': 'Retention Time'},
            {'field': 'Formula', 'title': 'Formula', 'headerTooltip': 'Chemical Formula'},
        ]

        # Create table
        self.table = pn.widgets.Tabulator(
            pagination=None,
            height=800,
            selectable='checkbox',
            header_align='center',
            layout='fit_columns',
            sizing_mode='stretch_width',
            show_index=False,
            theme='midnight',
            hidden_columns=['Compounds ID', 'abs_fc']
        )
        self.table.columns = self.table_columns

        # Plot pane
        self.plot_pane = pn.pane.Plotly(sizing_mode='stretch_width', min_height=900)

        # Set up callbacks
        self.comparison_select.param.watch(self._update_comparison, 'value')
        self.table.param.watch(self._update_selection, 'selection')
        self.search_input.param.watch(self._apply_filter, 'value')
        #self.select_all_button.on_click(self._select_all)
        self.clear_all_button.on_click(self._clear_all)
        self.reset_button.on_click(self._reset_app)
      
        # Initial update
        self._update_comparison(None)

    def _update_comparison(self, event):
        current_comp_idx = self.comparison_select.value
        self.comparison = current_comp_idx

        fc_col = self.comparisons[current_comp_idx].get("fold_change_col")
        pv_col = self.comparisons[current_comp_idx].get("p_value_col")

        # Prepare table data - include hidden columns for functionality
        table_data = self.df[[
            'Compounds ID', 'm/z', 'RT [min]', 'Formula', fc_col, pv_col,
        ]].dropna(subset=[fc_col]).copy()

        # Rename the columns in the dataframe
        table_data = table_data.rename(columns={
            'RT [min]': 'RT',
            pv_col: 'P-val',
            fc_col: 'Log2 Fold'
        })

        # Update table columns with the renamed columns
        self.table_columns = [
            {'field': 'm/z', 'title': 'm/z', 'headerTooltip': 'm/z value'},
            {'field': 'RT [min]', 'title': 'RT [min]', 'headerTooltip': 'Retention Time'},
            {'field': 'Formula', 'title': 'Formula', 'headerTooltip': 'Chemical Formula'},
            {'field': 'P-val', 'title': 'P-val', 'headerTooltip': 'P-value'},
            {'field': 'Log2 Fold', 'title': 'Log2 Fold', 'headerTooltip': 'Log2 Fold Change'},
        ]
        self.table.columns = self.table_columns

        # Sort by fold change magnitude
        table_data['abs_fc'] = abs(table_data['Log2 Fold'])
        table_data = table_data.sort_values('abs_fc', ascending=False)

        # Store the full dataset for filtering
        self.df_filtered = table_data

        # Reset search filter
        self.search_input.value = ''
        self.table.value = table_data
        self.table.selection = []

        # Update plot - need to use original column names for the plot
        self.plot_pane.object = generate_plot(self.df, current_comp_idx, self.comparisons)

    def _apply_filter(self, event):
        query = event.new.strip()
        if not query:
            self.table.value = self.df_filtered
            return

        # Filter only the 'm/z' column
        mask = self.df_filtered['m/z'].astype(str).str.contains(query)
        self.table.value = self.df_filtered[mask]

    def _update_selection(self, event):
        selected_ids = []
        if self.table.selection:
            selected_rows = [self.table.value.iloc[i] for i in self.table.selection]
            selected_ids = [str(row['Compounds ID']) for row in selected_rows]

        # Update plot with selection
        self.plot_pane.object = generate_plot(
            self.df, self.comparison, self.comparisons, selected_ids)

    #def _select_all(self, event):
    #    self.table.selection = list(range(len(self.table.value)))
  
    def _clear_all(self, event):
        self.table.selection = []
  
    def _reset_app(self, event):
        # Reset search filter
        self.search_input.value = ''
        # Reset comparison select to default
        self.comparison_select.value = 0
        # Clear table selection
        self.table.selection = []
        # Reload table data from the original dataframe and update plot
        self._update_comparison(None)
  
    def panel(self):

        # Button row
        button_row = pn.Row(
            #self.select_all_button,
            self.clear_all_button,
            self.reset_button,
            sizing_mode='stretch_width'
        )

        # Control panel with title, dropdown, search bar, buttons and table
        control_panel = pn.Column(
            "## Compound Explorer",
            self.comparison_select,
            self.search_input,
            button_row,
            self.table,
            styles={'background': '#606060'},
            css_classes=['custom-panel'],
            width=600,
            margin=0
        )

        # Main panel with plot
        main_panel = pn.Column(
            "## Volcano Plot",
            self.plot_pane,
            styles={'background': '#606060'},
            css_classes=['custom-panel'],
            min_width=550,
            margin=0,
            sizing_mode='stretch_width'
        )

        return pn.Row(
            control_panel,
            main_panel,
            styles={'background': '#606060'},
            sizing_mode='stretch_width',
            margin=0
        )

def notify_webhook():
        if ENV_CHECK != "DEV":
            webhook_url = "https://maker.ifttt.com/trigger/sidebar/json/with/key/" + IFTTT_KEY
            payload = {
                "value1": "Deploying Website"
            }
            headers = {
                "Content-Type": "application/json"
            }

            try:
                response = requests.post(webhook_url, json=payload, headers=headers)
                response.raise_for_status()
                #print("✅ IFTTT webhook triggered: Website deployed")
            except requests.exceptions.RequestException as e:
                print(f"[WARN] Webhook notify failed: {e}")
                if response is not None:
                    print(f"[DEBUG] Response status code: {response.status_code}")

# Main entry point
def main():
    app1 = VolcanoApp("ReSpleen.csv", "by_distance_named.csv", "ReSpleen.csv")
    app2 = VolcanoApp("ReKidney.csv", "ReKidney_by_distance_named.csv", "ReKidney.csv")
    app3 = VolcanoApp("ReLiver.csv", "ReLiver_by_distance_named.csv", "ReLiver.csv")

    tabs = pn.Tabs(
        ("Spleen", app1.panel()),
        ("Kidney", app2.panel()),
        ("Liver", app3.panel()),
    )

    if (ENV_CHECK == "DEV"):
        port = 4603
    else:
        port = 80
  
    pn.serve(tabs, port=port, websocket_origin=['*'])

if __name__ == "__main__":
    notify_webhook()
    main()