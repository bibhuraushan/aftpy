import os
import re
from urllib.parse import urlparse

def is_url(string):
    """
    Checks if the given string is a valid URL.

    Parameters
    ----------
    string : str
        The input string to check.

    Returns
    -------
    bool
        True if the string is a valid URL, False otherwise.
    """
    url_pattern = re.compile(
        r'^(https?|ftp)://[^\s/$.?#].[^\s]*$', re.IGNORECASE
    )
    return isinstance(string, str) and re.match(url_pattern, string) is not None


def is_path(string):
    """
    Checks if the given string is a valid file path.

    Parameters
    ----------
    string : str
        The input string to check.

    Returns
    -------
    bool
        True if the string is a valid file path, False otherwise.
    """
    return isinstance(string, str) and os.path.isfile(string)


def is_dir(string):
    """
    Checks if the given string is a valid directory path.

    Parameters
    ----------
    string : str
        The input string to check.

    Returns
    -------
    bool
        True if the string is a valid directory, False otherwise.
    """
    return isinstance(string, str) and os.path.isdir(string)


def is_email(string):
    """
    Checks if the given string is a valid email address.

    Parameters
    ----------
    string : str
        The input string to check.

    Returns
    -------
    bool
        True if the string is a valid email address, False otherwise.
    """
    email_pattern = re.compile(
        r'^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$'
    )
    return isinstance(string, str) and re.match(email_pattern, string) is not None


def is_int(string):
    """
    Checks if the given string represents an integer.

    Parameters
    ----------
    string : str
        The input string to check.

    Returns
    -------
    bool
        True if the string represents an integer, False otherwise.
    """
    return isinstance(string, str) and string.isdigit()


def is_float(string):
    """
    Checks if the given string represents a floating-point number.

    Parameters
    ----------
    string : str
        The input string to check.

    Returns
    -------
    bool
        True if the string represents a floating-point number, False otherwise.
    """
    float_pattern = re.compile(r'^-?\d+(\.\d+)?$')
    return isinstance(string, str) and re.match(float_pattern, string) is not None


def is_hex(string):
    """
    Checks if the given string represents a hexadecimal number.

    Parameters
    ----------
    string : str
        The input string to check.

    Returns
    -------
    bool
        True if the string represents a hexadecimal number, False otherwise.
    """
    hex_pattern = re.compile(r'^0x[0-9a-fA-F]+$')
    return isinstance(string, str) and re.match(hex_pattern, string) is not None


# Example Usage
if __name__ == "__main__":
    print("I am just a tool!")             # True