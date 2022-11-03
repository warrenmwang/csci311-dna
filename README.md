# CSCI311 DNA Project - Team 3
## Group Members: Nate Ahearn, Sam Baldwin, Nishant Shrestha, Warren Wang

### How To Use
Tested and written for Python `3.10.6` and above. Application is terminal/command line based, one can run in their bash interpreter. Additionally, in an interpreter with python runnability function, such as Intellij or VSCode, running from thhe IDE is also possible. <br>
`python ./main.py`

### Description
Running the code will present the user with a request for a file that contains the sequences that the query sequence will be compared to. If the file that is being selected is within the same directory as the .py file, only the filename is necessary, otherwise the full path to the file is necessary. 
The user will then be prompted with a request for the query file, containing the sequence that is to be compared. The same logic as above for providing the file name applies. 
If there is an issue in selecting either file, the program will continually request the names of files until valid files are given. 

Once the files have been provided, the program asks for input from the user, in the form of an integer from 0-3, to determine which of the algorithms shoudl be run on the two provided files. 
The algorithm will then print the comparison for each of the sequences within the sequences file against the query file, and the corresponding return value from the algorithm. 

After the algorithm that is chosen is finished running, the user is prompted if they want to run another algorithm or exit the program.
