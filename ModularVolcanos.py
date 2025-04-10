import pandas as pd
import plotly.graph_objects as go
import numpy as np
import json
import warnings
import re
import os


def clean_column_names(columns):
    return columns.str.strip().str.replace('"', '', regex=False).str.replace("'", '', regex=False)


def clean_cell_values(series):
    series = series.astype(object)
    series = series.where(series.notna(), '')
    series = series.astype(str)
    try:
        return series.replace(r'^\s*["\']{0,2}\s*$', np.nan, regex=True).str.strip(' "\'')
    except AttributeError:
        return series


def robust_load_csv(filepath, expected_columns=None):
    try:
        df = pd.read_csv(
            filepath,
            encoding="ISO-8859-1",
            engine="python",
            sep=",",
            quotechar='"',
            on_bad_lines="skip"
        )
    except Exception as e:
        raise RuntimeError(f"[ERROR] Could not load {filepath}: {e}")

    df.columns = clean_column_names(df.columns)
    df = df.loc[:, ~df.columns.duplicated()]

    for col in df.columns:
        df[col] = clean_cell_values(df[col])

    if expected_columns and not expected_columns.issubset(set(df.columns)):
        missing = expected_columns - set(df.columns)
        raise ValueError(f"[ERROR] Missing expected columns: {missing}")

    return df


def apply_fallback_names(df):
    for col in ['Name', 'Formula', 'Calc. MW']:
        if col not in df.columns:
            df[col] = np.nan
        else:
            df[col] = clean_cell_values(df[col])
    df['Name'] = df['Name'].fillna(df['Formula']).fillna(df['Calc. MW'])
    return df


def ensure_column(df, col, default=''):
    if col not in df.columns:
        df[col] = default
    else:
        df[col] = clean_cell_values(df[col])
    return df


def generate_volcano_plot(csv_file, mapping_file, distance_file, output_html):
    csv_key = os.path.basename(csv_file)

    # Load distance results
    distance_df = robust_load_csv(distance_file, expected_columns={'Compounds ID', 'Calc. MW', 'Name'})
    distance_df = apply_fallback_names(distance_df)
    gold_ids = distance_df['Compounds ID'].dropna().astype(str).tolist()

    # Load raw CSV data
    raw_data = robust_load_csv(csv_file)
    raw_data = apply_fallback_names(raw_data)
    raw_data = ensure_column(raw_data, 'm/z')
    raw_data = ensure_column(raw_data, 'RT [min]')
    raw_data['Compounds ID'] = raw_data['Compounds ID'].astype(str)
    raw_data['Gold'] = raw_data['Compounds ID'].isin(gold_ids)

    # Load column mapping
    with open(mapping_file) as f:
        mapping_data = json.load(f)

    comparisons = mapping_data.get(csv_key, [])
    if not comparisons:
        raise ValueError(f"No comparisons found under '{csv_key}' in {mapping_file}")

    # Build plot
    fig = go.Figure()
    dropdown_buttons = []

    for i, entry in enumerate(comparisons):
        fc_col = entry.get("fold_change_col")
        pv_col = entry.get("p_value_col")
        title = entry.get("title", f"Comparison {i+1}")

        if fc_col not in raw_data.columns or pv_col not in raw_data.columns:
            warnings.warn(f"Skipping comparison due to missing column(s): {fc_col}, {pv_col}")
            continue

        raw_data[fc_col] = pd.to_numeric(clean_cell_values(raw_data[fc_col]), errors='coerce')
        raw_data[pv_col] = pd.to_numeric(clean_cell_values(raw_data[pv_col]), errors='coerce')
        raw_data['-Log10(P-value)'] = -np.log10(raw_data[pv_col])

        sig_up = (raw_data[fc_col] > 0.5) & (raw_data[pv_col] < 0.05)
        sig_down = (raw_data[fc_col] < -0.5) & (raw_data[pv_col] < 0.05)

        color_map = np.where(raw_data['Gold'], 'gold',
                      np.where(sig_up, 'green',
                      np.where(sig_down, 'red', 'blue')))

        hover_text = (
            "Name:\t" + raw_data['Name'].astype(str) + "<br>" +
            "Formula:\t" + raw_data['Formula'].astype(str) + "<br>" +
            "Calc. MW:\t" + raw_data['Calc. MW'].astype(str) + "<br>" +
            "m/z:\t" + raw_data['m/z'].astype(str) + "<br>" +
            "RT [min]:\t" + raw_data['RT [min]'].astype(str)
        )

        fig.add_trace(go.Scatter(
            x=raw_data[fc_col],
            y=raw_data['-Log10(P-value)'],
            mode='markers',
            marker=dict(color=color_map, opacity=0.9),
            hovertext=hover_text,
            hoverinfo='text',
            visible=(i == 0),
            name=title,
            hoverlabel=dict(font_size=16, font_family="Arial", bgcolor="red", bordercolor="black")
        ))

        dropdown_buttons.append(dict(
            method="update",
            label=title,
            args=[{"visible": [j == i for j in range(len(comparisons))] + [True]*4},
                  {"title": title,
                   "xaxis": {"title": fc_col},
                   "yaxis": {"title": f"-Log10({pv_col})"}}]
        ))

    # Add legend traces
    legend_traces = [
        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=12, color='gold'),
                   showlegend=True, name='Top 50 Distance'),
        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=12, color='green'),
                   showlegend=True, name='Significant Upregulated'),
        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=12, color='red'),
                   showlegend=True, name='Significant Downregulated'),
        go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=12, color='blue'),
                   showlegend=True, name='Insignificant')
    ]

    for trace in legend_traces:
        fig.add_trace(trace)

    fig.update_layout(
        updatemenus=[dict(
            buttons=dropdown_buttons,
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=1.15,
            xanchor="right",
            y=1.15,
            yanchor="top"
        )],
        title=comparisons[0]["title"],
        xaxis_title=comparisons[0]["fold_change_col"],
        yaxis_title=f"-Log10({comparisons[0]['p_value_col']})",
        plot_bgcolor="lightslategray",
        paper_bgcolor="lightslategray",
        font=dict(color="black"),
        shapes=[
            dict(type="line", x0=-8, x1=8, y0=-np.log10(0.05), y1=-np.log10(0.05),
                 line=dict(color="Black", dash="dash")),
            dict(type="line", x0=-0.5, x1=-0.5, y0=0, y1=1, yref="paper",
                 line=dict(color="Black", dash="dash")),
            dict(type="line", x0=0.5, x1=0.5, y0=0, y1=1, yref="paper",
                 line=dict(color="Black", dash="dash"))
        ]
    )

    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)
    fig.write_html(output_html)
    print(f"✅ Volcano plot generated and saved as '{output_html}'")