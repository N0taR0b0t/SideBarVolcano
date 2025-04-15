# volcano_plot.py
import plotly.graph_objects as go
import numpy as np

def generate_plot(df, comparison_idx, comparisons, selected_ids=None):
    # Get comparison info
    entry = comparisons[comparison_idx]
    fc_col = entry.get("fold_change_col")
    pv_col = entry.get("p_value_col")
    title = entry.get("title", f"Comparison {comparison_idx+1}")

    # Create figure
    fig = go.Figure()

    # Determine opacity based on selection
    if selected_ids is None or len(selected_ids) == 0:
        opacity = np.ones(len(df))
    else:
        opacity = np.where(df['Compounds ID'].isin(selected_ids), 1.0, 0.005)

    # Determine colors
    color_map = np.where(df['Gold'], 'gold',
                  np.where(df[f'{fc_col}_sig_up'], 'green',
                  np.where(df[f'{fc_col}_sig_down'], 'red', 'blue')))

    # Create hover text
    hover_text = (
        "Name: " + df['Name'].astype(str) + "<br>" +
        "Formula: " + df['Formula'].astype(str) + "<br>" +
        "m/z: " + df['m/z'].astype(str) + "<br>" +
        "RT [min]: " + df['RT [min]'].astype(str) + "<br>" +
        "P-value: " + df[pv_col].apply(lambda x: f"{x:.2e}")
    )

    # Add main scatter plot
    fig.add_trace(go.Scatter(
        x=df[fc_col],
        y=df[f'-Log10({pv_col})'],
        mode='markers',
        marker=dict(color=color_map, opacity=opacity),
        hovertext=hover_text,
        hoverinfo='text',
        name=title,
        hoverlabel=dict(font_size=16, font_family="Arial", bgcolor="white", bordercolor="black")
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

    # Medium dark gray background with white text/grid
    fig.update_layout(
        title=title,
        xaxis_title=f"Log2 Fold Change: {fc_col}",
        yaxis_title=f"-Log10(P-value): {fc_col}",
        plot_bgcolor="#4f4f4f",
        paper_bgcolor="#4f4f4f",
        font=dict(color="white"),
        shapes=[
            dict(type="line", x0=-8, x1=8, y0=-np.log10(0.05), y1=-np.log10(0.05),
                line=dict(color="white", dash="dash")),
            dict(type="line", x0=-0.5, x1=-0.5, y0=0, y1=1, yref="paper",
                line=dict(color="white", dash="dash")),
            dict(type="line", x0=0.5, x1=0.5, y0=0, y1=1, yref="paper",
                line=dict(color="white", dash="dash"))
        ],
        legend=dict(font=dict(color="white"), bgcolor="#5f5f5f")
    )

    fig.update_xaxes(showgrid=True, gridcolor="white", gridwidth=0.5, color="white")
    fig.update_yaxes(showgrid=True, gridcolor="white", gridwidth=0.5, color="white")

    return fig