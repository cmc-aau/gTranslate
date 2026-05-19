import logging
import os
import sys
from functools import lru_cache

class __GTranslateCommonConfig:

    def __init__(self):
        # Internal settings used for logging.
        self.LOG_TASK = 21
        self._env_model_path = None

        self.logger = logging.getLogger('timestamp')


    @property
    def ENV_MODEL_PATH(self):
        if self._env_model_path is None:
            raw_path = os.getenv('GTRANSLATE_MODEL_PATH')

            if raw_path:
                self._env_model_path = os.path.expandvars(raw_path)
            else:
                self.logger.warning("The 'GTRANSLATE_MODEL_PATH' environment variable is not defined.")
        return self._env_model_path


    @property
    def ADA_MULTI_CLASS(self):
        return 'ada_multi_class.pkl.gz'

    @property
    def DT_MULTI_CLASS(self):
        return 'dt_multi_class.pkl.gz'

    @property
    def KNN_MULTI_CLASS(self):
        return 'knn_multi_class.pkl.gz'

    @property
    def XGB_MULTI_CLASS(self):
        return 'xgb_multi_class.pkl.gz'

    @property
    def MLP_MULTI_CLASS(self):
        return 'mlp_multi_class.pkl.gz'

    @property
    def LABEL_ENCODER(self):
        return 'label_encoder.pkl.gz'

# Export the class for import by other modules
@lru_cache(maxsize=1)
def __get_config():
    return __GTranslateCommonConfig()

CONFIG = __get_config()