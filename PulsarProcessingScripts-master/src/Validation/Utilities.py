"""
This code runs on python 2.4 or later.

By Rob Lyon <robert.lyon@cs.man.ac.uk>

+-----------------------------------------------------------------------------------------+
+                       PLEASE RECORD ANY MODIFICATIONS YOU MAKE BELOW                    +
+-----------------------------------------------------------------------------------------+
+ Revision |   Author    | Description                                       |    DATE    +
+-----------------------------------------------------------------------------------------+

 Revision:0    Rob Lyon    Initial version of the code.                        01/02/2014 

"""

# Python 2.4 imports.
import traceback
import sys

# ******************************
#
# CLASS DEFINITION
#
# ******************************

class Utilities(object):
    """
    Provides utility functions used when computing scores. Replaces the
    makeprofile.py file.
    
    """
    
    # ******************************
    #
    # Constructor.
    #
    # ******************************
    
    def __init__(self,debugFlag):
        self.debug = debugFlag
        
    # ******************************
    #
    # Functions.
    #
    # ******************************
    
    def appendToFile(self,path,text):
        """
        Appends the provided text to the file at the specified path.
        
        Parameters:
        path    -    the path to the file to append text to.
        text    -    the text to append to the file.
        
        Returns:
        N/A
        """
        
        destinationFile = open(path,'a')
        destinationFile.write(str(text))
        destinationFile.close()
    
    # ******************************
    
    def fileExists(self,path):
        """
        Checks a file exists, returns true if it does, else false.
        
        Parameters:
        path    -    the path to the file to look for.
        
        Returns:
        True if the file exists, else false.
        """
        
        try:
            fh = open(path)
            fh.close()
            return True
        except IOError:
            return False
    
    # ******************************************************************************************
            
    def format_exception(self,e):
        """
        Formats error messages.
        
        Parameters:
        e    -    the exception.
        
        Returns:
        
        The formatted exception string.
        """
        exception_list = traceback.format_stack()
        exception_list = exception_list[:-2]
        exception_list.extend(traceback.format_tb(sys.exc_info()[2]))
        exception_list.extend(traceback.format_exception_only(sys.exc_info()[0], sys.exc_info()[1]))
        
        exception_str = "Traceback (most recent call last):\n"
        exception_str += "".join(exception_list)
        
        # Removing the last \n
        exception_str = exception_str[:-1]
        
        return exception_str
    
    # ******************************************************************************************
            