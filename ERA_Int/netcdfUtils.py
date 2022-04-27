""" This module contains code for reading data from NetCDF files and
    intepreting metadata """
    
def getAttribute(var, attName, default=None):
    """ Gets the value of the given attribute as a string.  Returns the
    given default if the attribute is not defined. This is a useful
    "helper" method, which avoids AttributeErrors being raised if the
    attribute isn't defined.  If no default is specified, this function
    returns None if the attribute is not found. """
    if attName in var.ncattrs():
        return var.getncattr(attName)
    else:
        return default

def getTitle(var):
    """ Returns a title for a variable in the form "name (units)"
    The name is taken from the standard_name if it is provided, but
    if it is not provided, it is taken from the id of the variable
    (i.e. var._name).  If the units are not provided, the string
    "no units" is used instead. """
    standardName = getAttribute(var, 'standard_name', var._name)
    units        = getAttribute(var, 'units',         'no units')
    return "%s (%s)" % (standardName, units)

#######################################################################################
#####  The following functions test to see if coordinate variables represent geographic
#####  or time axes
#######################################################################################

def isLongitudeVar(coordVar):
    """ Given a coordinate variable (i.e. a NetCDF Variable object), this returns
    True if the variable holds values of longitude, False otherwise """
    # In the Climate and Forecast conventions, longitude variables are indicated
    # by their units (not by their names)
    units = getAttribute(coordVar, 'units')
    # There are many possible options for valid longitude units
    if units in ['degrees_east', 'degree_east', 'degree_E', 'degrees_E', 'degreeE', 'degreesE']:
        return True
    else:
        return False
    
def isLatitudeVar(coordVar):
    """ Given a coordinate variable (i.e. a NetCDF Variable object), this returns
    True if the variable holds values of latitude, False otherwise """
    # In the Climate and Forecast conventions, latitude variables are indicated
    # by their units (not by their names)
    units = getAttribute(coordVar, 'units')
    # There are many possible options for valid latitude units
    if units in ['degrees_north', 'degree_north', 'degree_N', 'degrees_N', 'degreeN', 'degreesN']:
        return True
    else:
        return False
    
def isVerticalVar(coordVar):
    """ Given a coordinate variable (i.e. a NetCDF Variable object), this returns
    True if the variable represents a vertical coordinate, False otherwise """
    # In the Climate and Forecast conventions, vertical coordinates are indicated
    # by units of pressure, or by the presence of a "positive" attribute.
    units = getAttribute(coordVar, "units")
    # First we look for units of pressure.  (There may be more possible pressure units
    # than are used here.)
    if units in ['Pa', 'hPa', 'pascal', 'Pascal']:
        return True
    else:
        # We don't have units of pressure, but perhaps we have a "positive" attribute
        positive = getAttribute(coordVar, 'positive')
        if positive in ['up', 'down']:
            return True
            
    # If we've got this far, we haven't satisfied either of the conditions for a
    # valid vertical axis
    return False
    
def isTimeVar(coordVar):
    """ Given a coordinate variable (i.e. a NetCDF Variable object), this returns
    True if the variable represents a time coordinate, False otherwise """
    # In the Climate and Forecast conventions, time coordinates are indicated
    # by units that conform to the pattern "X since Y", e.g. "days since 1970-1-1 0:0:0".
    # For simplicity, we just look for the word "since" in the units.  A complete
    # implementation should check this more thoroughly.
    units = getAttribute(coordVar, 'units')
    if units is None:
        # There are no units, so this can't be a time coordinate variable
        return False
    # The "find()" function on strings returns the index of the first match of the given
    # pattern.  If no match is found, find() returns -1.
    if units.find("since") >= 0:
        return True
    else:
        return False


#######################################################################################
#####  The following functions find geographic and time coordinate axes
#####  for data variables.
#######################################################################################

# As you can see, there is a lot of repetition in these functions - they all do basically
# the same thing.  There is a way to avoid this repetition, but it involves a technique
# that you may not be familiar with (i.e. passing functions as arguments to other functions)
# - see me if you want to know more.
    
def findLongitudeVar(nc, dataVar):
    """ Given a NetCDF Dataset object and a Variable object representing a data
    variable, this function finds and returns the Variable object representing the
    longitude axis for the data variable.  If no longitude axis is found, this returns
    None. """
    # First we iterate over the dimensions of the data variable
    for dim in dataVar.dimensions:
        # We get the coordinate variable that holds the values for this dimension
        coordVar = nc.variables[dim]
        # We test to see if this is a longitude variable, if so we return it
        if isLongitudeVar(coordVar):
            return coordVar
    # If we get this far we have not found the required coordinate variable
    return None

def findLatitudeVar(nc, dataVar):
    """ Given a NetCDF Dataset object and a Variable object representing a data
    variable, this function finds and returns the Variable object representing the
    latitude axis for the data variable.  If no latitude axis is found, this returns
    None. """
    for dim in dataVar.dimensions:
        coordVar = nc.variables[dim]
        if isLatitudeVar(coordVar):
            return coordVar
    return None

def findVerticalVar(nc, dataVar):
    """ Given a NetCDF Dataset object and a Variable object representing a data
    variable, this function finds and returns the Variable object representing the
    vertical axis for the data variable.  If no vertical axis is found, this returns
    None. """
    for dim in dataVar.dimensions:
        coordVar = nc.variables[dim]
        if isVerticalVar(coordVar):
            return coordVar
    return None

def findTimeVar(nc, dataVar):
    """ Given a NetCDF Dataset object and a Variable object representing a data
    variable, this function finds and returns the Variable object representing the
    time axis for the data variable.  If no time axis is found, this returns
    None. """
    for dim in dataVar.dimensions:
        coordVar = nc.variables[dim]
        if isTimeVar(coordVar):
            return coordVar
    return None

def isPositiveUp(zVar):
    """ Given a vertical coordinate variable, this function returns true if the
    values on the vertical axis increase upwards. For vertical axes based on pressure,
    the values increase downward, so this returns False.  If the axis is not based
    on pressure, the value of the "positive" attribute (which can be "up" or "down")
    is used instead. """
    units = getAttribute(zVar, "units")
    # First we look for units of pressure.  (If we find them, this is a pressure axis
    # and the values increase downward)
    if units in ['Pa', 'hPa', 'pascal', 'Pascal']:
        return False
    else:
        # We don't have units of pressure, but perhaps we have a "positive" attribute
        positive = getAttribute(zVar, 'positive')
        if positive == 'up':
            return True
        else:
            return False
			
			
##################################################################################
##### PROJECT FILES ############################################################

######### QUESTION 1 ########################################################
# just a copy of findNearestIndex from utils.py    
def findNearestIndex(vals, target):
    """ Searches through vals for the target value, returning the index
    of the value in vals that is closest numerically to the target value.
    If more than one value in vals is equally close to the target, the index
    of the first value will be returned. """
    minIndex = -1
    minDiff = None
    # Loop over all the values in the given list
    for i in range(len(vals)):
        # Find the absolute difference between this value and the target
        diff = abs(vals[i] - target)
        # If this is the first time, or if the difference is smaller than
        # the smallest difference found so far, remember the difference and
        # the index
        if minDiff is None or diff < minDiff:
            minDiff = diff
            minIndex = i
    return minIndex

# now for Latitude
def findNearestLatIndex(nc,dataVar,LatVal):
    ''' This function first finds the latitude coordinate variable with 
    latitude values in given data variable object (e.g. temperature). Then 
    it searches for nearest index to the given latitude value and returns 
    the index. If more than one value in values is equally close to the LatVal, 
    the index of the first value will be returned. 
    If there is no latitude variable a ValueError is raised. 
    Note that it is assumed that Latitude values are always between -90 and 90.'''
    # First find the coordinate variable representing latitude in dataVar
    # ps. is variable object!
    coordVar=findLatitudeVar(nc, dataVar)
    
    if coordVar is None:
        # check if no latitude values were found
        raise ValueError("Cannot search for latitude value if latitude axis is not present")
    else:
        # if latitude values were found
        # now find the values of coordVar
        vals = coordVar[:] # Read all lat values
        # now find the nearest index
        index=findNearestIndex(vals, LatVal)
        return index

# for Vertical axis
def findNearestZIndex(nc,dataVar,ZVal):
    ''' This function first finds the vertical coordinate variable with 
    height/depth/pressure values in given data variable object (e.g. temperature). 
    Then it searches for nearest index to the given vertical value and returns 
    the index. If more than one value in values is equally close to the ZVal, 
    the index of the first value will be returned. 
    If there is no vertical variable a ValueError is raised.'''
    # First find the coordinate variable representing vertical axis in dataVar
    # ps. is variable object!
    coordVar=findVerticalVar(nc, dataVar)
    
    if coordVar is None:
        # check if no Z values were found
        raise ValueError("Cannot search for latitude value if latitude axis is not present")
    else:
        # if Z values were found
        # now find the values of coordVar
        vals = coordVar[:] # Read all Z values
        # now find the nearest index
        index=findNearestIndex(vals, ZVal)
        return index
		
# for Longitude axis

# define a function that can convert between longitude values 
# from 0 to 360 and longitude values from -180 to 180
def lon180(lon):
    '''this function converts a number that is larger than 180 to a negative number 
    to obtain longitude values between -180 and 180 degrees, where -180=180.'''
    lon = lon % 360 # -25%360=335 !!!!! 40%360=40 !!!!
    if lon > 180:
        return lon - 360
    else:
        return lon  
	
# then find nearest index
def findNearestLonIndex(nc,dataVar,LonVal):
    ''' This function first finds the longitude coordinate variable with 
    longitude values in given data variable object (e.g. temperature). Then 
    it searches for nearest index to the given longitude value and returns 
    the index. If more than one value in values is equally close to the LonVal, 
    the index of the first value will be returned. 
    If there is no longitude variable a ValueError is raised. 
    Note that longitudes can be either between -180 and 180 degrees or 
    between 0 and 360 degrees, e.g. if LonVal is -25 we seek for either 
    -25 or 335 degrees.'''
    # First find the coordinate variable representing longitude in dataVar
    # ps. is variable object!
    coordVar=findLongitudeVar(nc, dataVar)
    
    if coordVar is None:
        # check if no longitude values were found
        raise ValueError("Cannot search for longitude value if longitude axis is not present")
    else:
        # if longitude values were found
        # now find the values of coordVar
        vals = coordVar[:] # Read all lon values
        
        # convert value of longitude to [-180,180] whatever the value was before
        LonVal=lon180(LonVal) 
        # now find appropriate index!
        if min(vals)<0:
		    # we have regime when lon is between -180 and 180
            # now find the nearest index
            index=findNearestIndex(vals, LonVal)
            return index
        else: 
            # min(vals) >= 0
            LonVal1 = LonVal % 360 # -25%360=335 !!!!! 40%360=40 !!!!
            # so now we are between 0 and 360 degrees
            # calculate nearest index
            index=findNearestIndex(vals, LonVal1)
            return index

#################### QUESTION 2.1 ###########################################

def getTitleZ(nc,var, direction,latOrlon):
    ''' Returns a title for a variable in the form "EW/NS section 
    of name (units) at certain degrees latitude/longitude" to 1 
    decimal place accurate. The name is taken from the standard_name 
    if it is provided, but if it is not provided, it is taken from 
    the id of the variable (i.e. var._name).  If the units are not 
    provided, the string "no units" is used instead. Also, if, 
    for example, EW section lacks the latitude axis it writes 
    "at unknown latitude". '''
    # note: var= variable object
    standardName = getAttribute(var, 'standard_name', var._name)
    units        = getAttribute(var, 'units',         'no units')
    if direction=="EW":
        # we have a section at certain latitude
        if findLatitudeVar(nc, var)==None:
            return "EW section of %s (%s) at unknown latitude" % (standardName, units)
        else:
            # we found a latitude that corresponds to given value
            return "EW section of %s (%s) at %.1f degrees latitude" % (standardName, units,latOrlon)
    else:
        # direction=="NS" we have a section at certain longitude
        if findLongitudeVar(nc, var)==None:
            return "NS section of %s (%s) at unknown longitude" % (standardName, units)
        else:
            # we found a longitude that corresponds to given value
            return "NS section of %s (%s) at %.1f degrees longitude" % (standardName, units,latOrlon)

#################### QUESTION 2.2 ###########################################

def getName(var,coordinate,value):
    ''' Returns units, value & name for a variable in the form "value units name"
    The name is taken from the standard_name if it is provided, but
    if it is not provided, it is taken from the id of the variable
    (i.e. var._name). If the units are not provided, the string "no units" 
    is used instead. If variable object does not exist 
    (its 'value' is None) then it returns "unknown coordinate". 
    Done for coordinate variables.'''
    if var is None:
        return "unknown %s" % (coordinate)
    else:
        standardName = getAttribute(var, 'standard_name', var._name)
        units        = getAttribute(var, 'units',         'no units')
        return "%.1f %s %s" % (value,units,standardName)

def getTitleT(nc,dataVar, lat,lon,Z):
    ''' Returns a title for a variable in the form "Timeseries 
    of name (units) at certain degrees latitude, certain degrees 
    longitude and/or certain elevation" to 1 decimal place accurate. 
    The name is taken from the standard_name 
    if it is provided, but if it is not provided, it is taken from 
    the id of the variable (i.e. var._name).  If the units are not 
    provided, the string "no units" is used instead. Also, if, 
    for example, the file lacks the latitude axis it writes 
    "at unknown latitude". '''
    # note: dataVar= variable object
    standardName1 = getAttribute(dataVar, 'standard_name', dataVar._name)
    units1        = getAttribute(dataVar, 'units',         'no units')
    # find variable objects (if they exist)
    ZVar=findVerticalVar(nc, dataVar) # find variable object for Z
    latVar=findLatitudeVar(nc, dataVar) # find variable object for lat
    lonVar=findLongitudeVar(nc, dataVar) # find variable object for lon
    lonName=getName(lonVar,"longitude",lon) # find the name and units of lon variable
    zName=getName(ZVar,"elevation",Z) # find the name and units of Z variable
    latName=getName(latVar,"latitude",lat) # find the name and units of lat variable

    return "Timeseries of %s (%s) at %s, %s and %s" % (standardName1, units1,latName,lonName,zName)
	










