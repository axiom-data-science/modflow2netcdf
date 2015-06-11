"""

Provides functionality to validate different types of multi-dimensional data and log validation
results to a validation log file.

"""

import os
import numpy as np


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class ValidationLog(object):
    """
    Data validation log class.  Logs various errors, warnings, and informational messages in an easy to read
    tabbed validation log file.

    Parameters
    ----------
    log_file_path : string
        Log file path including the name of the log file.  Both relative and full paths are accepted.

    Methods
    -------
        write_data_type_not_supported
        write_test_header
        write_test_data_header
        write_test_data_without_test
        write_test_data_missing
        write_data_element_failure
        write_test_failure
        close
    See Also
    --------

    Notes
    -----

    Examples
    --------

    >>> import data_validation_engine
    >>> dve_log = data_validation_engine.ValidationLog(log_file_path)

    """
    def __init__(self, log_file_path):
        self.logfile = open(log_file_path, 'w')

    def write_data_type_not_supported(self, data_type):
        self.logfile.write('\tWARNING: Unable to test data of type %s.  Data type not currently supported\n' % data_type)
        self.logfile.flush()

    def write_test_header(self, test_name):
        self.logfile.write('\n*** Testing project: %s ***\n' % test_name)
        self.logfile.flush()

    def write_test_result(self, result, test_name):
        if result :
            self.logfile.write('\t\tTest %s succeeded.\n' % test_name)
        else:
            self.logfile.write('\t\tTest %s failed.\n' % test_name)

    def write_test_data_header(self, test_data_name):
        self.logfile.write('\tTesting data: %s\n' % test_data_name)
        self.logfile.flush()

    def write_test_data_without_test(self, test_data_name):
        self.logfile.write('\tWARNING: No test data found to test: %s\n' % test_data_name)
        self.logfile.flush()

    def write_test_data_missing(self, test_data_name):
        self.logfile.write('\tERROR: Expected test data not found: %s\n' % test_data_name)
        self.logfile.flush()

    def write_data_element_failure(self, file_name, element_number, data_expected, data_to_validate):
        self.logfile.write('\t\tERROR: Data verification failed for data element comparing against file %s.  '
                           'Element number %s did not match.\n' % (file_name, element_number))
        if is_float(data_to_validate):
            self.logfile.write('\t\t\tFound:    %f\n' % data_to_validate)
        else:
            self.logfile.write('\t\t\tFound:    %s\n' % str(data_to_validate))

        if is_float(data_to_validate):
            self.logfile.write('\t\t\tExpected: %f\n' % data_expected)
        else:
            self.logfile.write('\t\t\tFound:    %s\n' % str(data_to_validate))

        self.logfile.flush()

    def write_test_failure(self, test_name, error_message=None):
        if error_message is not None:
            self.logfile.write('ERROR: Internal test %s failed with error message: %s\n.' % (test_name, error_message))
        else:
            self.logfile.write('ERROR: Internal test %s failed with an unexpected error.\n' % test_name)
        self.logfile.flush()

    def close(self):
        self.logfile.close()


class IntegerWildcardString(object):
    """
    Represents a string that contains multiple wildcard characters.  Each wildcard character must be an integer.
    Wildcard characters are in the format of:
        <WildCardChar>X
    where X is a single character that represents a unique wildcard that
    has an integer value.

    Parameters
    ----------
    wild_card_string : string
        Representation of a string using wild card characters
    wild_card_char : string
        Special character that indicates the next character is a wildcard.

    Methods
    -------
        resolve_string
    See Also
    --------

    Notes
    -----

    Examples
    --------

    >>> import data_validation_engine
    >>> wildcardchar = data_validation_engine.IntegerWildcardString('This string has #1 two types of wildcard #2 characters #1.', '#')
    >>> wildcardchar.resolve_string()

    """

    def __init__(self, wild_card_string, wild_card_char):
        self.string_parts = []
        self.wild_cards = []
        self.wild_card_vals = []

        raw_string_parts = wild_card_string.split(wild_card_char)
        for part in raw_string_parts:
            if len(self.string_parts) == 0:
                self.string_parts.append(part)
            else:
                self.string_parts.append(part[1:])
                self.wild_cards.append(part[0])

    # Resolves wildcard values specified in dtParams from a string.
    # Input -   input_string: String to be resolved
    # Output -  Returns True if the string can be resolved, false otherwise.
    def resolve_string(self, input_string):
        """
        Determines if the passed in string is an instance of the wild card string represented by the instance of this \
        class.  If the string is an instance of this wild card string, resolve_string determines all the values of
        wild cards in the string and stores them in self.wild_card_vals.

        Parameters
        ----------
        input_string : String
            String to be represented as this wild card string.

        Returns
        -------
        boolean : True if the input_string is an instance of this wild card string, false otherwise

        Examples
        --------

        """

        self.wild_card_vals = []

        # Locate all strings
        match_string = input_string
        index = -1
        for string_part in self.string_parts:
            if index != -1:
                match_string = self._resolve_wildcard(match_string, string_part)
            if match_string is None or not match_string.startswith(string_part):
                return False
            else:
                match_string = match_string[len(string_part):]

            index += 1

        return len(match_string) == 0

    def _resolve_wildcard(self, match_string, string_part):
        # Resolve largest wildcard string so that the next match string matches
        match_found_index = -1
        wc_break_index = 0
        wild_card_part = match_string[:wc_break_index]
        while is_int(wild_card_part) and wc_break_index + 1 < len(match_string):
            if match_string[wc_break_index+1:].startswith(string_part):
                match_found_index = wc_break_index
            wc_break_index += 1
            wild_card_part = match_string[:wc_break_index]
        # Check for valid wild card
        if match_found_index == -1:
            return None
        # Update with new wild card
        self.wild_card_vals.append(int(match_string[:match_found_index]))
        return match_string[match_found_index+1:]


class TestProjectData(object):
    def __init__(self,strDataName):
        self.strDataName = strDataName
        self.arrData = []


# Configuration information for a specific data set
class TestProjectDataType(object):
    def __init__(self,strDataName,strDataType,strDataFilePath,strDataFileName,strDataDelimiter):
        self.strDataName = strDataName
        self.strDataType = strDataType
        self.strDataFilePath = strDataFilePath
        self.strDataFileName = strDataFileName
        self.strDataDelimiter = strDataDelimiter
        self.bHasBeenTested = False
        self.arrData = []


def no_conversion(data):
    return data


# Parent class for all data tests
class DataValidationEngine(object):
    def __init__(self, data_type_supported, validation_folder_path, log_file, data_delimiter, data_error_threshold=0.0):
        self.file_names = []
        self.resolved_params = []

        # Set up converter for supported data types
        if data_type_supported == 'float':
            self.data_convert = float
        elif data_type_supported == 'integer':
            self.data_convert = int
        else:
            self.data_convert = no_conversion

        self.data_type_supported = data_type_supported
        self.validation_folder_path = validation_folder_path
        self.log_file = log_file
        self.data_delimiter = data_delimiter
        self.data_error_threshold = data_error_threshold

    def validate(self, file_name, data_to_validate, params=None):
        """
        Method of abstract base class.  Subclasses validate a data set against an expected data set.

        Parameters
        ----------
        file_name - String
            Data file name to look for.  Supports wild card characters.
        data_to_validate - Numpy array
            validation data stored in a numpy array of any dimension
        params -
            any additional parameters needed to validate a given data type

        Returns
        -------
        boolean : True if the data is valid, false otherwise

        Examples
        --------

        """

        # Abstract method should not be called
        raise Exception('DataValidationEngines abstract method validate called.')

    # Resolves a file name with wildcards
    # Input
    #       wildcard_file_name - Name of the file, with wildcard characters, to search for
    # Output
    #       Returns the file name of the resolved file and stores the parameters for that
    #       file in dtResolvedParams.  Returns None if not found.
    def resolve_names(self, wildcard_file_name, wildcard_char):
        resolved_params = []

        # Break wildcard string into array
        intwildcard = IntegerWildcardString(wildcard_file_name, wildcard_char)

        # Loop through files in data folder
        for file_name in os.listdir(self.validation_folder_path):
            # If file name indicates this is a data file to be processed
            if intwildcard.resolve_string(file_name):
                resolved_params.append((intwildcard.wild_card_vals, file_name))

        return resolved_params

    def _data_file_compare(self, data_file_name, validation_data_iterator):
        success = True
        element_number = 1
        data_file = open(data_file_name, 'r')
        for line in data_file:
            if self.data_delimiter == '':
                line_array = line.split()
            else:
                line_array = line.split(self.data_delimiter)

            for data_expected in line_array:
                data_expected_casted = self.data_convert(data_expected)
                # Either data must be exactly as expected or, if an error threshold is set,
                # must be within the error threshold
                data_to_validate = validation_data_iterator.next()[1]
                if (self.data_error_threshold != 0.0 and abs(data_to_validate - data_expected_casted) >
                    self.data_error_threshold) or (self.data_error_threshold == 0.0 and
                    data_to_validate != data_expected_casted):
                    # Report error to log file
                    self.log_file.write_data_element_failure(data_file_name, element_number, data_expected_casted, data_to_validate)
                    success = False
                element_number += 1
        return success


class ValidateDataArray(DataValidationEngine):
    @classmethod
    def is_registrar_for(cls, data_type):
        return data_type == 'DataArray'

    def validate(self, file_name, data_to_validate, params=None):
        validation_data_iterator = np.ndenumerate(data_to_validate)
        file_path = os.path.join(self.validation_folder_path, file_name)
        return self._data_file_compare(file_path, validation_data_iterator)


class ValidateTimeDependentDataArray(DataValidationEngine):
    @classmethod
    def is_registrar_for(cls, data_type):
        return data_type == 'TimeDependentDataArray'

    # Subclassed from DataValidationEngine to validate time dependent data.  See parent class description.
    def validate(self, file_name, data_to_validate, wildcard_char=None):
        success = True
        # Get data shape
        validation_data_iterator = np.ndenumerate(data_to_validate)
        # Find all data files that fit the data file name
        resolved_params = self.resolve_names(file_name, wildcard_char)
        if len(resolved_params) > 0:
            # Sort in order of integer variables found
            resolved_params = sorted(self.resolved_params, key=lambda params: (params[0]))
            # Loop through each file to be tested and test with next part of data set
            for params in self.resolved_params:
                file_path = os.path.join(self.validation_folder_path, params[1])
                if not self._data_file_compare(file_path, validation_data_iterator):
                    success = False
        return success


# Class factory for all test engines.  Returns the correct test engine for a given data type.
def test_engine_factory(data_type, data_check, validation_folder_path, log_file, data_delimiter, data_error_threshold=0.0):
    for cls in DataValidationEngine.__subclasses__():
        if cls.is_registrar_for(data_check):
            return cls(data_type, validation_folder_path, log_file, data_delimiter, data_error_threshold)
    return None
