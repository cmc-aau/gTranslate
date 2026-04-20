import os
import logging
import pandas as pd
import plotly.express as px

from gtranslate.biolib_lite.common import make_sure_path_exists


class FeaturePlotter:
    """Generates an interactive HTML 2D/3D scatter plot from the classifier features."""

    def __init__(self, tsv_path: str, out_file: str, highlight_file: str = None):
        self.tsv_path = tsv_path
        self.html_path = out_file
        self.highlight_file = highlight_file
        # make sure the self.html_path ends with .html
        if not self.html_path.endswith('.html'):
            self.html_path += '.html'
        self.logger = logging.getLogger('timestamp')

    def generate_html(self):
        """Reads the TSV and generates the standalone HTML interactive plot."""
        if not os.path.exists(self.tsv_path):
            self.logger.warning(f"Cannot generate plot: {self.tsv_path} not found.")
            return

        try:
            # 1. Load the data
            df = pd.read_csv(self.tsv_path, sep='\t')

            if df.empty:
                self.logger.warning("Feature file is empty. Skipping plot generation.")
                return

            # 2. Process the highlight file (if provided)
            highlight_genomes = set()
            if self.highlight_file and os.path.exists(self.highlight_file):
                with open(self.highlight_file, 'r') as f:
                    highlight_genomes = {line.strip() for line in f if line.strip()}

            if highlight_genomes:
                # Create a column that flags if the genome is in the highlight list
                df['Highlight'] = df['user_genome'].apply(
                    lambda x: 'Highlighted' if x in highlight_genomes else 'Normal'
                )
                self.logger.info(f"Highlighting {len(highlight_genomes)} genomes in the plot.")

            # 3. Setup Colors and Categories
            color_col = None
            categories = []
            cat_orders = None
            color_map = {}  # Used to force 'Highlighted' to be gold

            if 'predicted_tln_table' in df.columns:
                df['predicted_tln_table'] = df['predicted_tln_table'].astype(str)
                # Create a dedicated column for plotting so we don't overwrite the actual TT data
                df['plot_color_group'] = df['predicted_tln_table']
                color_col = 'plot_color_group'
                categories = ['11', '4', '25']
            elif highlight_genomes:
                # Fallback if there is no TT column but highlights exist
                df['plot_color_group'] = 'Normal'
                color_col = 'plot_color_group'
                categories = ['Normal']
                color_map['Normal'] = '#636efa'  # Default Plotly blue

            if highlight_genomes:
                # Override the color group for highlighted genomes
                df.loc[df['Highlight'] == 'Highlighted', color_col] = 'Highlighted'
                if 'Highlighted' not in categories:
                    categories.append('Highlighted')
                # Assign gold to the highlighted category
                color_map['Highlighted'] = 'gold'

            if color_col:
                cat_orders = {color_col: categories}

            # 4. Get the feature columns
            ignore_cols = {'user_genome', 'predicted_tln_table', 'Highlight', 'plot_color_group'}
            features = [col for col in df.columns if col not in ignore_cols]

            # Setup hover data to preserve the original TT info for highlighted dots
            hover_cols = []
            if 'predicted_tln_table' in df.columns:
                hover_cols.append('predicted_tln_table')

            if len(features) < 2:
                self.logger.warning("Not enough features to generate a plot (need at least 2).")
                return

            def create_buttons(dimension, is_3d=True):
                buttons = []
                for f in features:
                    if color_col:
                        trace_data = [df[df[color_col] == cat][f].values for cat in categories]
                    else:
                        trace_data = [df[f].values]

                    layout_update = {}
                    if is_3d:
                        layout_update[f'scene.{dimension}axis.title.text'] = f
                    else:
                        layout_update[f'{dimension}axis.title.text'] = f

                    # 'update' method accepts [data_updates, layout_updates]
                    buttons.append(dict(
                        method='update',
                        label=f,
                        args=[{dimension: trace_data}, layout_update]
                    ))
                return buttons

            # 5. Create the marker size slider
            size_steps = []
            for s in [2, 4, 6, 8, 10, 14, 20]:
                size_steps.append(
                    dict(
                        method='restyle',
                        label=str(s),
                        args=[{'marker.size': s}]
                    )
                )

            size_slider = [dict(
                active=2,
                currentvalue={"prefix": "Marker Size: "},
                pad={"t": 50},
                steps=size_steps
            )]

            if len(features) >= 3:
                fig = px.scatter_3d(
                    df, x=features[0], y=features[1], z=features[2],
                    color=color_col, hover_name='user_genome',
                    hover_data=hover_cols,  # Added to keep original TT visible on hover
                    category_orders=cat_orders, title="gTranslate Feature Explorer",
                    color_discrete_map=color_map  # Forces the gold color mapping
                )

                fig.update_layout(
                    margin=dict(l=250),
                    sliders=size_slider,
                    updatemenus=[
                        dict(buttons=create_buttons('x', True), direction="down", x=-0.25, y=0.9, xanchor="left",
                             yanchor="top", showactive=True),
                        dict(buttons=create_buttons('y', True), direction="down", x=-0.25, y=0.7, xanchor="left",
                             yanchor="top", showactive=True),
                        dict(buttons=create_buttons('z', True), direction="down", x=-0.25, y=0.5, xanchor="left",
                             yanchor="top", showactive=True)
                    ],
                    annotations=[
                        dict(text="<b>X-Axis</b>", x=-0.25, y=0.92, xref="paper", yref="paper", showarrow=False,
                             xanchor="left", yanchor="bottom"),
                        dict(text="<b>Y-Axis</b>", x=-0.25, y=0.72, xref="paper", yref="paper", showarrow=False,
                             xanchor="left", yanchor="bottom"),
                        dict(text="<b>Z-Axis</b>", x=-0.25, y=0.52, xref="paper", yref="paper", showarrow=False,
                             xanchor="left", yanchor="bottom")
                    ]
                )

            else:
                fig = px.scatter(
                    df, x=features[0], y=features[1],
                    color=color_col, hover_name='user_genome',
                    hover_data=hover_cols,  # Added to keep original TT visible on hover
                    category_orders=cat_orders, title="gTranslate Feature Explorer",
                    color_discrete_map=color_map  # Forces the gold color mapping
                )

                fig.update_layout(
                    margin=dict(l=250),
                    sliders=size_slider,
                    updatemenus=[
                        dict(buttons=create_buttons('x', False), direction="down", x=-0.25, y=0.9, xanchor="left",
                             yanchor="top", showactive=True),
                        dict(buttons=create_buttons('y', False), direction="down", x=-0.25, y=0.7, xanchor="left",
                             yanchor="top", showactive=True)
                    ],
                    annotations=[
                        dict(text="<b>X-Axis</b>", x=-0.25, y=0.92, xref="paper", yref="paper", showarrow=False,
                             xanchor="left", yanchor="bottom"),
                        dict(text="<b>Y-Axis</b>", x=-0.25, y=0.72, xref="paper", yref="paper", showarrow=False,
                             xanchor="left", yanchor="bottom")
                    ]
                )

            make_sure_path_exists(os.path.dirname(self.html_path))
            fig.write_html(self.html_path, include_plotlyjs="cdn")
            self.logger.info(f"Interactive HTML dashboard generated: {self.html_path}")

        except Exception as e:
            self.logger.error(f"Failed to generate feature explorer plot: {e}")