###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import logging
import sys
import traceback

from gtranslate import __author__, __copyright__, __version__
from gtranslate.biolib_lite.exceptions import BioLibError
from gtranslate.biolib_lite.logger import logger_setup
from gtranslate.cli import get_main_parser
from gtranslate.exceptions import GTranslateExit, GTranslateException
from gtranslate.main import OptionsParser


def print_help():
    print('''\

              ...::: gTranslate v%s :::...

  Tools:
    detect_table  -> Detect the genetic translation table (GTT) used in prokaryotic organisms.
    generate_plot -> Generate an interactive HTML dashboard to explore the features used for GTT prediction.

  Testing:
    test          -> Validate detection of the genetic translation table (GTT) used in prokaryotic organisms.
    check_install -> Check the installation of gTranslate.
    
  Training:
    ground_truth   -> Determine the ground truth for genomes based on their taxonomic classification.
    build_features -> Generate feature vectors for training models.
    fit_models     -> Train models based on the feature vector tables.
    training_wf    -> Run the full training workflow in one command ( ground_truth -> build_features -> fit_models ).

  Use: gtranslate <command> -h for command specific help 
    ''' % __version__)


def main():
    # -------------------------------------------------
    # get and check options
    args = None
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f"gtranslate: version {__version__} {__copyright__} {__author__}")

        sys.exit(0)
    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    else:
        args = get_main_parser().parse_args()

    # setup logger
    logger_setup(args.out_dir if hasattr(args, 'out_dir') and args.out_dir else None,
                 "gtranslate.log", "gTranslate", __version__, False,
                 hasattr(args, 'debug') and args.debug)
    logger = logging.getLogger('timestamp')

    # -------------------------------------------------
    # do what we came here to do
    try:
        gt_parser = OptionsParser(__version__,
                                  args.out_dir if hasattr(args, 'out_dir') and args.out_dir else None)
        gt_parser.parse_options(args)
    except SystemExit:
        logger.error('Controlled exit resulting from early termination.')
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error('Controlled exit resulting from interrupt signal.')
        sys.exit(1)
    except GTranslateExit as e:
        if len(str(e)) > 0:
            logger.error('{}'.format(e))
        logger.error('Controlled exit resulting from an unrecoverable error or warning.')
        sys.exit(1)
    except (GTranslateException, BioLibError) as e:
        msg = 'Controlled exit resulting from an unrecoverable error or warning.\n\n'
        msg += '=' * 80 + '\n'
        msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
        msg += '  MESSAGE: {}\n'.format(e)
        msg += '_' * 80 + '\n\n'
        msg += traceback.format_exc()
        msg += '=' * 80
        logger.error(msg)
        sys.exit(1)
    except Exception as e:
        msg = 'Uncontrolled exit resulting from an unexpected error.\n\n'
        msg += '=' * 80 + '\n'
        msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
        msg += '  MESSAGE: {}\n'.format(e)
        msg += '_' * 80 + '\n\n'
        msg += traceback.format_exc()
        msg += '=' * 80
        logger.error(msg)
        sys.exit(1)


if __name__ == '__main__':
    main()
