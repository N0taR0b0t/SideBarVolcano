# Root Cause of the Memory Leak
The memory leak is caused by repeated creation and assignment of large Plotly Figure objects to the self.plot_pane.object attribute within each VolcanoApp, without proper cleanup or reuse.

In VolcanoApp.__init__(), _update_comparison() is called once per app instance (3 total).
Each call to _update_comparison() calls generate_plot(...), which creates a brand-new Plotly Figure.
This figure is then assigned to self.plot_pane.object.
Panel's Plotly pane under the hood wraps the Plotly figure in a Bokeh model. When .object is reassigned:

The old figure may not be immediately garbage-collected, especially if:
Panel keeps internal references (which it often does to preserve UI state),
The Plotly figure includes large numpy arrays, color maps, or binary buffers,
Or if the assignment happens frequently (e.g., in a callback or re-rendering loop).
Over time, if new Figure objects keep replacing old ones without explicitly dereferencing the prior ones or invoking gc.collect(), Python's garbage collector delays cleanup, especially if circular references or C-level memory (e.g. Plotly's JSON serialization buffers) are involved.

Even though _update_comparison() is only called once at startup for each app:

Three full Plotly figures are instantiated and retained simultaneously.
If the VolcanoApp is tabbed via pn.Tabs, each app and its corresponding figure remains in memory, even if not visible.
Panel's live server (via pn.serve) does not automatically prune hidden tabs or reclaim their memory.
Plotly figures are large — they embed full hover text, color arrays, shape objects, etc.
