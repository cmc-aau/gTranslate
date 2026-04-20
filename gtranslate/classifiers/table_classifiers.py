import logging
import joblib
import pandas as pd

from gtranslate.config.common import CONFIG
from gtranslate.exceptions import GTranslateExit


class GenericTableClassifier:
    """A generic table classifier that uses a list of column classifiers to classify tables."""
    def __init__(self, scaler_path :str, classifier_path : str, tt_mapping : dict):
        self.scaler_path = scaler_path
        self.classifier_path = classifier_path
        self.tt_mapping = tt_mapping
        self.logger = logging.getLogger('timestamp')

    def load_model(self, model_path):
        """Load a model from the given path, with error handling."""
        try:
            return joblib.load(model_path)
        except FileNotFoundError:
            self.logger.error("Model file not found: %s", model_path)
            raise GTranslateExit(f"Model file not found: {model_path}")
        except Exception as e:
            self.logger.error("Error loading model from %s: %s", model_path, e)
            raise GTranslateExit(f"Error loading model from {model_path}: {e}")


    def scale_data(self, data : pd.DataFrame) -> pd.DataFrame:
        """Scale the input DataFrame to match the trained scaler."""
        scaler = self.load_model(self.scaler_path)
        # make sure temp_df has the same columns as the scaler
        temp_df = data.reindex(columns=scaler.feature_names_in_, fill_value=0)
        # transform the genome info
        temp_df_scaled = scaler.transform(temp_df)
        # Convert the scaled array back into a DataFrame with the original column names
        temp_df_scaled = pd.DataFrame(temp_df_scaled, columns=temp_df.columns, index=temp_df.index)
        return temp_df_scaled

    def classify_table(self, data : pd.DataFrame) -> tuple[int, float]:
        """Classify the input data, returning (class, probability)."""
        classifier = self.load_model(self.classifier_path)
        # I want to make the classifier non verbose
        if hasattr(classifier, 'verbose'):
            classifier.verbose = 0
        # The classifier is train on a specific subset of columns so we need to make sure the data has the same columns
        trimmed_data = data.reindex(columns=classifier.feature_names_in_, fill_value=0)
        probabilities = classifier.predict_proba(trimmed_data)
        predictions = classifier.predict(trimmed_data)
        return predictions[0], max(probabilities[0])

    def scale_and_classify(self, data : pd.DataFrame):
        """Scale and classify the data."""
        scaled_data = self.scale_data(data)
        return self.classify_table(scaled_data)

    def predict_translation_table(self, data: pd.DataFrame):
        """Get the translation table for the list of genomes."""
        binary_value_prediction, probability = self.scale_and_classify(data)
        return self.tt_mapping[binary_value_prediction], probability

class Classifier_4_11(GenericTableClassifier):
    """A classifier that distinguishes between translation tables 4 and 11."""

    def __init__(self,classifier_path :str=None, scaler_path : str=None):
        super().__init__(scaler_path or CONFIG.SCALER_4_11,
                         classifier_path or CONFIG.CLASSIFIER_4_11,
                         {0:4, 1:11})


class Classifier_25(GenericTableClassifier):
    """A classifier that distinguishes between translation table 25 and 4."""

    # we can add the path to the classifiers and scalers None by default
    def __init__(self, classifier_path :str=None, scaler_path : str=None):
        super().__init__(scaler_path or CONFIG.SCALER_25,
                         classifier_path or CONFIG.CLASSIFIER_25,
                         {1:4, 0:25})



