import pandas as pd
import bs4 as bs
import numpy as np
import requests
import datetime as dt
import pathlib
import os
import drms
from sunpy.time import parse_time

home = pathlib.Path.home()


def wso(start_date="1976-01-01", end_date=dt.datetime.today(), outfile=None):
    """
    Fetch and process the polar magnetic field data from the Wilcox Solar Observatory (WSO).

    The function retrieves the polar field data, processes it, and returns a pandas DataFrame.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with columns:
        - "Time": Timestamps (ISO 8601 format)
        - "PolarN": Northern polar field strength
        - "PolarS": Southern polar field strength
    """
    # Fetch the webpage content
    url = "http://wso.stanford.edu/Polar.html"
    response = requests.get(url)
    soup = bs.BeautifulSoup(response.content, 'html5lib')

    # Initialize lists for storing time and polar field data
    time, pn, ps = [], [], []

    # Save raw data to a local file
    with open(home / "temp_wso_polarfield_raw.temp", "w") as fl:
        for line in soup.pre.text.splitlines()[1:]:
            fl.write(line + "\n")

    # Process the raw data file
    with open(home / "temp_wso_polarfield_raw.temp", "r") as fl:
        for line in fl:
            # Parse each line
            _array = line.split()
            if len(_array) == 0:
                continue  # Skip empty lines

            try:
                # Parse time in ISO format
                time.append(dt.datetime.strptime(_array[0], "%Y:%m:%d_%Hh:%Mm:%Ss").isoformat())

                # Parse polar field strengths
                pn.append(float(_array[1].removesuffix("N")))
                ps.append(float(_array[2].removesuffix("S")))
            except (ValueError, IndexError):
                # Handle parsing errors or missing data
                pn.append(np.nan)
                ps.append(np.nan)

    # Create a DataFrame from the parsed data
    df = pd.DataFrame({"Time": pd.to_datetime(time), "PolarN": pn, "PolarS": ps}).set_index("Time")
    start_date=pd.to_datetime(parse_time(start_date).to_datetime())
    end_date=pd.to_datetime(parse_time(end_date).to_datetime())
    df = df[(df.index >= start_date) & (df.index <= end_date)]
    os.remove(home / "temp_wso_polarfield_raw.temp")

    return df


def hmi(start_date=None, end_date=None, outfile=None, verbose=True):
    """
    Retrieve HMI mean magnetic field data between specified dates using JSOC DRMS.

    Parameters
    ----------
    start_date : str or datetime.datetime, optional
        Start date for data retrieval (default: 1 year before `end_date`).
    end_date : str or datetime.datetime, optional
        End date for data retrieval (default: today's date).
    outfile : str, optional
        Path to save the output DataFrame as a CSV file (if provided).
    verbose : bool, optional
        If True, prints query progress and results (default is True).

    Returns
    -------
    pandas.DataFrame
        DataFrame containing HMI mean magnetic field data with datetime index.

    Raises
    ------
    ValueError
        If the date inputs are not valid datetime or ISO format strings.
    """
    # Handle end_date
    end_date = parse_time(end_date).to_datetime() if end_date else dt.datetime.today()

    # Handle start_date
    start_date = parse_time(start_date).to_datetime() if start_date else end_date - dt.timedelta(days=365)

    # Ensure valid date range
    if start_date >= end_date:
        raise ValueError("`start_date` must be earlier than `end_date`.")

    # Date formatting
    fmt = "%Y.%m.%d_%H:%M:%S_TAI"

    # Generate yearly time ranges
    time_range = pd.date_range(start_date, end_date, freq="YS")  # Yearly Start

    # Ensure at least two time points for iteration
    if len(time_range) == 1 or time_range[-1] < end_date:
        time_range = time_range.append(pd.DatetimeIndex([end_date]))

    # Initialize DRMS client once
    c = drms.Client()

    # Query in yearly chunks
    df_list = []
    for t1, t2 in zip(time_range[:-1], time_range[1:]):
        if verbose:
            print(f"Querying {t1.date()} to {t2.date()}...", end=" ")

        t_start = t1.strftime(fmt)
        t_end = t2.strftime(fmt)

        # DRMS query
        qstr = f'hmi.meanpf_720s[{t_start}-{t_end}][? QUALITY=0 ?]'
        keys = 'T_REC,CRLT_OBS,CAPN1,CAPS1,CAPN2,CAPS2'
        _df = c.query(qstr, key=keys)

        if _df.empty:
            if verbose:
                print("No data found.")
            continue

        if verbose:
            print(f"Retrieved {_df.shape[0]} records.")

        df_list.append(_df)

    # Concatenate results
    if not df_list:
        raise ValueError("No data retrieved for the specified date range.")

    df = pd.concat(df_list).drop_duplicates()

    # Set datetime index
    df.set_index('T_REC', inplace=True)
    df.index = pd.to_datetime(df.index, format=fmt)
    df.index.name = "Time"

    # Optionally save to file
    if outfile:
        df.to_csv(outfile)
        if verbose:
            print(f"Data saved to {outfile}")

    return df
