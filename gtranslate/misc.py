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
import os

from gtranslate.biolib_lite.execute import check_dependencies
from gtranslate.biolib_lite.logger import colour
from gtranslate.exceptions import GTranslateExit


class Misc(object):

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger('timestamp')


    def checkfile(self, file_path, file_name):
        """Check that a file exists, output the result to the logger.

        Returns
        -------
        bool
            True if the folder exists, False otherwise.
        """
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            self.logger.info("Check file {} ({}): {}".format(
                file_name, file_path, colour('OK', ['bright'], fg='green')))
            return True
        else:
            self.logger.warning("Check file {} ({}): {}".format(
                file_name, file_path, colour('MISSING', ['bright'], fg='red')))
            return False

    def checkfolder(self, folder_path, folder_name):
        """Check that a folder exists, output the result to the logger.

        Returns
        -------
        bool
            True if the folder exists, False otherwise.
        """
        if os.path.isdir(folder_path) and len(os.listdir(folder_path)) > 0:
            self.logger.info("Check folder {} ({}): {}".format(
                folder_name, folder_path, colour('OK', ['bright'], fg='green')))
            return True
        else:
            self.logger.warning("Check folder {} ({}): {}".format(
                folder_name, folder_path, colour('MISSING', ['bright'], fg='red')))
            return False



    def check_install(self):
        """Check that all reference files exist.

        Returns
        -------
        bool
            True if the installation is complete, False otherwise.
        """

        # Check that all programs are on the system path.
        self.logger.info(f'Checking that all third-party software are on the system path:')
        names = {'prodigal'}
        for name in sorted(names):
            on_path = False
            try:
                on_path = on_path or check_dependencies([name], exit_on_fail=False)
            except:
                pass
            if on_path:
                self.logger.info("         |-- {:16} {}".format(
                    name, colour('OK', ['bright'], fg='green')))
            else:
                self.logger.info("         |-- {:16} {}".format(
                    name, colour('NOT FOUND', ['bright'], fg='yellow')))

        # Assume this was successful unless otherwise observed.
        ok = True

        if not ok:
            raise GTranslateExit('Unexpected files were seen, or the reference package is corrupt.')
