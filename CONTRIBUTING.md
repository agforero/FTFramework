### How to helpfully contibute to the framework:
In order to reduce margin of error with the compilation tests, it would be helpful to have as many test subjects as possible. Thus, if you have working FORTRAN code ready in a directory, fork this repository, use `./addNew <directory>` on the directory, and send a pull request! 

Be careful not to submit any code with race conditions present. To check, add it to your fork, and use `./verify.py` from `/main/`.

### Another way to help:
As of now, I have `compare.py` set to detect GNU Make errors whenever a string fulfills any of the following conditions in `errorDeclaration()`:

```Python3
def errorDeclaration(line):                     # is the line make telling us that there's been an error
    # this function was made by referring to https://www.gnu.org/software/make/manual/html_node/Error-Messages.html
    banana = line.split()
    if len(banana) < 2: return False            # we need at least a make: and something after
    if banana[0] == "#": banana = banana[1:]    # shave off the pound sign if it's there
    if banana[0][:5] == "make:":                # according to link above, this will always be true unless the Makefile itself is botched
        if banana[1] == "***": return True      # should always be enough with my flags
        elif banana[-2] == "Error":             # if not, see if "Error" is in the line
            try: 
                int(banana[-1])                 # see if you can convert string after "Error" into an integer
                return True 
            except: pass                        # suppress the error message, even though this case is weird
        elif "Segfault" in banana: return True  # maybe it segfaulted? we could go on with elifs at this point
    return False
```

If you notice GNU Make outputting any error that does NOT fit the requirements above, please edit this function to include it, and submit a pull request. According to the documentation referenced in the first comment after the `def`, an error occuring in the compiler itself should always contain `make: ***`; but there have been examples seen where this isn't the case, at least in earlier GNU Make versions. You could just add more `elif` statements to account for your newfound error(s).
