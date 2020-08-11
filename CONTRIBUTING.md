### How to helpfully contibute to the framework:
In order to reduce margin of error with the compilation tests, it would be helpful to have as many test subjects as possible. Thus, if you have working FORTRAN code ready in a directory, fork this repository, use `./addNew <directory>` on the directory, and send a pull request! 

Be careful not to submit any code with race conditions present. To check, add it to your fork, and use `./verify.py` from `/main/`.
