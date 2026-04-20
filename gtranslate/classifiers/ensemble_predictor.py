import os
from collections import Counter

import joblib
import numpy as np

from gtranslate.config.common import CONFIG

class TTPredictor:
    def __init__(self,custom_model_path=None):
        """
        Initializes the predictor by loading models directly via the centralized CONFIG.
        """


        # Helper function to decide which path to use
        def resolve_path(default_config_path):
            if custom_model_path is not None:
                # Extract the filename (e.g., 'ada_multi_class.pkl.gz') and join it with the custom dir
                filename = os.path.basename(default_config_path)
                return os.path.join(custom_model_path, filename)
            return default_config_path

        # Load models using the resolved paths
        self.ada = joblib.load(resolve_path(CONFIG.ADA_MULTI_CLASS))
        self.dt = joblib.load(resolve_path(CONFIG.DT_MULTI_CLASS))
        self.knn = joblib.load(resolve_path(CONFIG.KNN_MULTI_CLASS))
        self.xgb = joblib.load(resolve_path(CONFIG.XGB_MULTI_CLASS))
        self.mlp = joblib.load(resolve_path(CONFIG.MLP_MULTI_CLASS))

        self.models = [self.ada, self.dt, self.knn, self.xgb, self.mlp]

        # Load label encoder
        self.label_encoder = joblib.load(resolve_path(CONFIG.LABEL_ENCODER))

    def predict_translation_table(self, df,custom_model_path=None):
        """
        Predicts using all 5 models, takes the majority vote, and decodes the label.
        """
        expected_features = self.ada.feature_names_in_
        warnings = []

        try:
            df_aligned = df[expected_features]
        except KeyError as e:
            raise ValueError(f"Your DataFrame is missing expected features: {e}")

        all_predictions = np.column_stack([model.predict(df_aligned) for model in self.models])
        decoded_preds = self.label_encoder.inverse_transform(all_predictions[0])
        decoded_preds_str = [str(p) for p in decoded_preds]

        vote_preds = Counter(decoded_preds)
        top_preds = vote_preds.most_common()

        best_class, max_votes = top_preds[0]

        if max_votes == 2:
            warnings.append(f"Low confidence: maximum model agreement was only {max_votes}/5.")

            recoding_votes = vote_preds.get('4', 0) + vote_preds.get('25', 0)

            if recoding_votes >= 3:
                # if 2x(TT4) vs 2x(TT25) Default to 4 because it is more biologically common, but add a warning
                # Pick the recoding table that had the most votes
                if vote_preds.get('4', 0) > vote_preds.get('25', 0):
                    best_class = '4'
                    warnings.append(f"Tie broken: Ensemble vote detects a recoding event with predictions {','.join(decoded_preds)}")
                elif vote_preds.get('4', 0) == vote_preds.get('25', 0):
                    best_class = '4'
                    warnings.append("Tie broken: 2 votes for TT4 and 2 votes for TT25. Defaulted to TT4.")
                else:
                    warnings.append(f"Tie broken: Ensemble vote detects a recoding event with predictions {','.join(decoded_preds)}")
                    best_class = '25'


        confidence_score = vote_preds.get(best_class, 0) / len(self.models)
        bonus_confidence = 0.0

        # Add 0.05 partial confidence for each vote going to the other recoding event
        if best_class == '4':
            bonus_confidence = vote_preds.get('25', 0) * 0.05
        elif best_class == '25':
            bonus_confidence = vote_preds.get('4', 0) * 0.05

        confidence_score += bonus_confidence

        ensemble_preds = {}


        for model, pred in zip(self.models, decoded_preds_str):
            # Extract the model's name
            if hasattr(model, 'steps'):
                model_name = list(model.named_steps.values())[-1].__class__.__name__
            else:
                model_name = type(model).__name__


            # Save the prediction to the dictionary, casting back to an integer for cleanliness
            # remove the 'Classifier' prefix from the model name for readability
            ensemble_preds[model_name] = int(pred)

        formatted_preds = {key.replace('Classifier', '').lower() + '_pred': value for key, value in
                       ensemble_preds.items()}

        feature_vector = df_aligned.iloc[0].to_dict()


        return best_class, confidence_score, warnings, formatted_preds,feature_vector
