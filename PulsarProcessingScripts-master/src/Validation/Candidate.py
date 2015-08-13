"""
Represents a pulsar candidate file (phcx file).

By Rob Lyon <robert.lyon@cs.man.ac.uk>

+-----------------------------------------------------------------------------------------+
+                       PLEASE RECORD ANY MODIFICATIONS YOU MAKE BELOW                    +
+-----------------------------------------------------------------------------------------+
+ Revision |   Author    | Description                                       |    DATE    +
+-----------------------------------------------------------------------------------------+

 Revision:0    Rob Lyon    Initial version of the code.                        03/02/2014 
 
 
"""

# ******************************
#
# CLASS DEFINITION
#
# ******************************

class Candidate:
    """
    Represents a pulsar candidate. This class is used as part of the validation
    code for the score generator application.
    
    """
    
    # ******************************
    #
    # INIT FUNCTION
    #
    # ******************************
    def __init__(self,name="Unknown",path=""):
        """
        Represents an individual Pulsar candidate.
        Expects an input parameter name which is the primary name for 
        this individual candidates - file name typically used. Also expects
        the path to the candidate file.
        
        """
        self.candidateName = name # Name of the candidate file, minus the path.
        self.candidatePath = path # The full path to the candidate.
        self.scores = []          # Stores all 22 candidate scores.
        self.label = "Unknown"    # The label this candidate has received, either "POSITIVE" or "NEGATIVE".
        
        # Some candidates may be 'special' in that they possess extreme
        # values for one of their scores. This could be a minimum or max
        # value. This variable stores the index of that special score
        # (if it exists) with respect to the scores[] array.
        self.specialScore=-1      
        
        # If this candidate does have a special score, this variable
        # is used to provide a single word description of why it is
        # 'special'. So for example if specialScore=1 then special="MAX"
        # would indicate this candidate has the highest score 1 known.
        self.special="None"       
        
    # ******************************
    #
    # UTILITY FUNCTIONS.
    #
    # ******************************
    
    def addScores(self,lineFromFile):
        """
        Adds the scores read in from the candidate .dat file to this object.
        
        Parameters:
        lineFromFile    -    the string text from the file.
        
        Returns:
        
        N/A
        
        """
        substrings = lineFromFile.split(",")

        counter = 1
        for s in substrings:
            if(s != "" or len(s)!=0):
                counter+=1
                self.scores.append(float(s))
            
    # ******************************
    # 
    # ******************************
    
    def getScore(self,index):
        """
        Obtains the specified score for this candidate. Compensates
        for zero indexing. So if score 1 is desired simply call
        getScore(1).
        
        Parameters:
        index    -    the index of the score to obtain.
        
        Returns:
        
        The floating point value of the desired score.
        """
        return float(self.scores[index-1])
    
    def getName(self):
        """
        Obtains the name of the candidate file, not the full path.
        
        
        Returns:
        
        The name of the candidate file.
        """
        return self.candidateName
    
    def getPath(self):
        """
        Obtains the full path to the candidate.
        
        
        Returns:
        
        The full path to the candidate.
        """
        return self.candidatePath
    
    def setLabel(self,l):
        """
        Sets the label describing this candidate, i.e. positive or negative.
        To be clear the input should either be l="POSITIVE" or l="NEGATIVE".
        
        Parameters:
        l    -    the label for this candidate, i.e. l="POSITIVE" or l="NEGATIVE".
        
        Returns:
        
        N/A
        """
        self.label = l
        
    def getLabel(self):
        """
        Gets the label describing this candidate, i.e. "POSITIVE" or "NEGATIVE".
        
        Parameters:
        N/A
        
        Returns:
        
        The string label describing this candidate.
        """
        return self.label
         
    def isPulsar(self):
        """
        Checks the label on this candidates, and determines if it
        represents a pulsar or not.
        
        Parameters:
        N/A
        
        Returns:
        
        True if this candidate represents a genuine pulsar, else False.
        """
        if(self.label=="POSITIVE"):
            return True
        else:
            return False
        
    def setSpecialScore(self,special):
        """
        Sets the value of the score which makes this candidate unusual,
        i.e. score 1 may be the maximum observed or the minimum observed.
        
        Parameters:
        special    -    the score which makes this candidate is unique.
        
        Returns:
        
        N/A
        """
        try:
            self.specialScore=int(special)
        except Exception as e: # catch *all* exceptions
            self.specialScore=-1
            
    def getSpecialScore(self):
        """
        Gets the value of the score which makes this candidate unusual.
        
        Parameters:
        
        N/A
        
        Returns:
        
        The integer value of the special score for this candidate.
        """
        try:
            return int(self.specialScore)
        except:
            return -1
    
    def setScores(self, data):
        """
        Sets the value of the scores for this candidate, stores them
        as an array of floating point values.
        
        Parameters:
        
        Data    -    the 22 candidate scores.
        
        Returns:
        
        N/A
        """
        
        self.scores=[float(i) for i in data]
        
    def setSpecial(self, s):
        """
        Sets the value of the special description. This should
        be either MAX or MIN. This would indicate along with the
        specialScore why this candidate is unusual, e.g.
        
        If specialScore = 5 and special= MAX then this candidate
        would be unusual since it has the maximum value for score
        5. Since we also have access to the candidate's true label,
        we could go further and say that it has the MAX score 5 value
        for the positive or the negative class (here positive means
        legitimate pulsar, negative RFI etc).
        
        Parameters:
        
        s    -    the string special description.
        
        Returns:
        
        N/A
        """
        self.special = str(s)
    
    def getSpecial(self):
        """
        Gets the value of the special description. This should
        be either MAX or MIN. This would indicate along with
        specialScore why this candidate is unusual, e.g.
        
        If specialScore = 5 and special= MAX then this candidate
        would be unusual since it has the maximum value for score
        5. Since we also have access to the candidate's true label,
        we could go further and say that it has the MAX score 5 value
        for the positive or the negative class (here positive means
        legitimate pulsar, negative RFI etc).
        
        Parameters:
        
        N/A
        
        Returns:
        
        Gets the string value of the special description.
        """
        return str(self.special)
        
    # ******************************
    # 
    # ******************************
    
    def __str__(self):
        """
        Overridden method that provides a neater string representation
        of this class. This is useful when writing these objects to a file
        or the terminal.
        
        """
            
        return self.candidateName + "," + self.candidatePath
    