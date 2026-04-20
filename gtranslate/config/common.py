import os
from functools import lru_cache

class __GTranslateCommonConfig:

    # Internal settings used for logging.
    LOG_TASK = 21

    def _get_model_path(self, filename):
        """Internal helper to construct absolute paths for classifier models."""
        current_dir = os.path.dirname(__file__)
        path = os.path.join(current_dir, '..', 'classifiers', 'classifier_models', filename)
        return os.path.abspath(path)

    @property
    def ADA_MULTI_CLASS(self):
        return self._get_model_path('ada_multi_class.pkl.gz')

    @property
    def DT_MULTI_CLASS(self):
        return self._get_model_path('dt_multi_class.pkl.gz')

    @property
    def KNN_MULTI_CLASS(self):
        return self._get_model_path('knn_multi_class.pkl.gz')

    @property
    def XGB_MULTI_CLASS(self):
        return self._get_model_path('xgb_multi_class.pkl.gz')

    @property
    def MLP_MULTI_CLASS(self):
        return self._get_model_path('mlp_multi_class.pkl.gz')

    @property
    def LABEL_ENCODER(self):
        return self._get_model_path('label_encoder.pkl.gz')

# Export the class for import by other modules
@lru_cache(maxsize=1)
def __get_config():
    return __GTranslateCommonConfig()

CONFIG = __get_config()