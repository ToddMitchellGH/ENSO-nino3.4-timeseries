def arclength( lat1, lon1, lat2, lon2, radius=None ):
    """Arc-length distance in km

    Assumes angles in degrees ( not radians ). 

    Todd Mitchell, April 2019"""

    if radius is None:
        radius = 6.37e3  # in km

    import numpy as np
    meanlat = np.mean( ( lat1, lat2 ) )
    rcosine   = radius * np.cos( np.deg2rad( meanlat ) )
    a = rcosine * ( lon2 - lon1 ) / 360  * 2 * np.pi
    b = radius  * ( lat2  - lat1 )   / 360 * 2 * np.pi

    return np.sqrt( a * a + b * b )

def fill_year( array, nperyr=None ):
    '''Append missing values ("NaN"s) to a time series / data file
    to make complete years of data.

    Input: 
    array      a time series or data matrix where time varies by row
    nperyer    (optional) number of records in a year; default 12

    Example: tpm.fill_year( ts )

     Todd Mitchell, June 2018 '''
    import numpy as np
    if nperyr is None:
        nperyr = 12
    nperyr *= 1.0
    if len( array.shape )==1: 
        nt = array.shape
        nt = nt[0]
        nt2 = np.ceil( nt/nperyr ) * nperyr
        nt2 = int( nt2 )
        temp = np.ones( ( nt2-nt ), dtype="int" ) * np.nan
        array = np.hstack( ( array, temp ) )
    else: 
        nt, nx = array.shape
        nt2 = np.ceil( nt/nperyr ) * nperyr
        nt2 = int( nt2 )
        temp = np.ones( ( nt2-nt, nx ), dtype="int" ) * np.nan
        array = np.vstack( ( array, temp ) )
    return( array )
def find_latlon( xgrid, ygrid, lat, lon ):
    '''Identify the x and y gridpoints that are closest to 
    the input lat and lon.  

    The code was written for Huancayo Peru, which is at a longitude of -75.21 
    and the input grid spans 275 to 330.  Be careful that the answer you 
    get makes sense ! 

    Huancayo -12.0668 latitude and -75.2103 longitude

    The algorithm assumes that xgrid, ygrid are vectors.

    The algorithm flags if 2 longitude or 2 latitude points are returnd.

    Todd Mitchell, January 2019'''

    import numpy as np

    temp = abs( ygrid - lat )
    yval = np.where( temp == min(temp) )
    if len(yval) > 1:
        print( len(yval), ' gridpoints nearest to ', lat )

    if lon<0:
        lon = lon + 360

    temp = abs( xgrid - lon )
    xval = np.where( temp == min(temp) )
    if len(xval) > 1: 
        print( len(xval), ' gridpoints nearest to ', lon )

    return{ 'yval': yval, 'xval': xval }
def space_longitudes():
    '''Adjust the aspect ratio of maps to take into account the 
    convergence of map meridians.
     Todd Mitchell, October 2018 '''
    import matplotlib.pyplot as plt
    import numpy as np
    xl, xr = plt.xlim()   # yields the x limits
    yb, yt = plt.ylim()
    plt.axes().set_aspect( aspect=( xr - xl ) / ( ( yt - yb ) * np.cos( np.radians( np.mean( [ yt, yb ] ) ) ) ) )
    return()
def threetotwo( array ):
    '''Turns a 3-dimensional array( nx, ny, nt ) into a 
    2-dimensional array ( nt, nx*ny ).

    Todd Mitchell, November 2018.'''
    import numpy as np
    nt, ny, nx = array.shape
    if min( [ nt, ny, nx ] )==1:
        print( 'This is a 2-dimensional array.' )
    array = np.transpose( array )
    array = np.reshape( array, ( nx*ny, nt ) )
    array = array.T
    return( array )
def write_ts( ts, yr1, yr2, yrfst=None ):
    '''write_ts( ts, yr1, yr2 ) writes a monthly timeseries in table form to stdio.

    This is the beginning of converting write_ts.m to python.

     Todd Mitchell, February 2019 '''

    import numpy as np
    import sys
    
    if yrfst is None:
        yrfst = yr1

    nyr = yr2 - yr1 + 1
    
    print( 'ts.shape yields', ts.shape )

    a = np.reshape( np.round( ts*10 ), ( nyr, 12 ) )
    a[ np.isnan(a) ] = -999
    b = np.arange( yr1, yr2+1 )
    b = np.expand_dims( b, axis=1 )
    b = np.concatenate( ( b, a ), axis=1 ).astype(int)
    np.savetxt( sys.stdout, b, fmt='%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d' )
def yearsmonths( yr1, yr2=None ):
    '''yearsmonths(yr1, yr2 ) returns a dictionary of monthly values of years,
    and calendar month for the input years

     yr2 optional, default is yr2 = yr1

     Output ã dictionary of monthly values of 
     years           year
     months          calendar month

     Example: a = tpm.yearsmonthsdays( 1970 ), a[ 'years' ] = 1970, 1970, ...

     Todd Mitchell, January 2019 '''

    import numpy as np
    if yr2 is None:
        yr2 = yr1

    nyr = yr2 - yr1 + 1

    years = np.arange( yr1, yr2+1 )
    years = np.expand_dims( years, axis=1 )
    ones = np.ones( ( 1, 12 ) )
    years = years * ones
    years = np.reshape( years, ( int(nyr)*12, 1 ) )

    months = np.arange( 1, 12+1 )
    months = np.expand_dims( months, axis=1 )
    ones = np.ones( ( 1, int(nyr) ) )
    months = ( months * ones ).T
    months = np.reshape( months, ( int(nyr)*12, 1 ) )

    return { 'years': years, 'months':months }
def yearsmonthsdays( yr1, yr2=None ):
    '''yearsmonthsdays(yr1, yr2 ) returns a dictionary of daily values of years,
    calendar month, day of month, and Julian day for the input years

     yr2 optional, default is yr2 = yr1

     Output ã dictionary of daily values of 
     years           year
     months          calendar month
     days            day of month
     jdays           Julian Day ( 1 - 365 ) or ( 1 - 366 for Leap Years )

     Example: a = tpm.yearsmonthsdays( 1970 ), a[ 'years' ] = 1970, 1970, ...

     Todd Mitchell, September 2017 '''
    import numpy as np
    if yr2 is None:
        yr2 = yr1
    ndays = np.array( [    # leap and 3 non-leap years, number of days in each calendar month, 
        [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ],
        [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ], 
        [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ], 
        [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ] ] )
    for year in np.arange( yr1, yr2+1 ):
        a = abs( np.mod( year-2000,4 ) )+1;   # a = 1 for a leap year
        for month in np.arange( 12 ):
            yearstemp = np.ones( ( ndays[a-1,month],1 ), dtype="int" ) * year
            monthstemp = np.ones( ( ndays[a-1,month],1 ), dtype="int" ) * (month+1)
            daystemp = np.arange( ndays[a-1,month-1] ).reshape( ( -1, 1 ) ) + 1
#           print( "year, month", year, month )
            if month==0 and year==yr1:
#               print( "Inside the month==0/year==yr1 branch." )
                years = yearstemp
                months = monthstemp
                days = daystemp
#               print( "Inside the month==0/year==yr1 branch.  days.shape", days.shape )
            else:
#               print( "Inside the else branch.  days.shape, daystemp.shape", days.shape, daystemp.shape )
                years = np.vstack( ( years, yearstemp ) )
                months = np.vstack( ( months, monthstemp ) )
                days = np.vstack( ( days, daystemp ) )
#               print( "Inside the else branch.  days.shape", days.shape )
#       print( "After the if/else code.  days.shape", days.shape )
        jdaystemp = np.arange(366)+1
        if a>1: jdaystemp = np.arange(365)+1
        jdaystemp = jdaystemp.reshape( ( -1, 1 ) )
        if year==yr1:
            jdays = jdaystemp
        else:
            jdays = np.vstack( ( jdays, jdaystemp ) )
        days = days.flatten(1)
        jdays = jdays.flatten(1)
        days = days.reshape( (-1, 1 ) )
        jdays = jdays.reshape( (-1, 1 ) )
#    return ( years, months, days, jdays )
    return { 'years': years, 'months':months, 'days':days, 'jdays':jdays }
