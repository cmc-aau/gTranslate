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


class GTranslateException(Exception):
    """ Base exception for all gTranslate exceptions thrown in this project. """

    def __init__(self, message=''):
        Exception.__init__(self, message)


class GTranslateExit(Exception):
    """Raised when gTranslate is to quietly exit."""

    def __init__(self, message=''):
        Exception.__init__(self, message)


class GenomeNameInvalid(GTranslateException):
    """ Thrown when a genome name contains characters which are not supported. """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)


class GenomeBatchfileMalformed(GTranslateException):
    """ Thrown when the format of the genome batchfile is malformed. """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)


class NoGenomesFound(GTranslateException):
    """ Thrown when no input genomes are found in a directory. """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)


class ReferenceFileMalformed(GTranslateException):
    """ Thrown when a reference file is malformed. """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)


class GenomeMarkerSetUnknown(GTranslateException):
    """ Thrown when the genome marker set is unknown (i.e. not ar53, or bac120). """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)


class InconsistentGenomeBatch(GTranslateException):
    """ Thrown when number of genomes in the identify directory is different than the number of genomes to process. """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)


class FileNotFound(GTranslateException):
    """ Thrown when a file is not found. """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)


class DirNotFound(GTranslateException):
    """ Thrown when a directory is not found. """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)


class ProdigalException(GTranslateException):
    """ Thrown when Prodigal returns a non-zero exit code. """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)

class GTranslateTestFailure(GTranslateException):
    """ Thrown when the gTranslate user test suite fails. """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)

class GTranslateArgsParsingConflict(GTranslateException):
    """ Thrown when the arguments are conflicting and or missing.  """

    def __init__(self, message=''):
        GTranslateException.__init__(self, message)
