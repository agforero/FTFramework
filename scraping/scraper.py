#!/usr/bin/env python3

import requests, sys, os
from bs4 import BeautifulSoup

def fix(name):
    if len(name) == 0: return ""
    elif name[-1] == "/": return ""
    else: return str(fix(name[:-1])) + str(name[-1])

def websiteName(url):
    try:
        for i in range(len(url)):
            if url[i:i+7] == "http://": # getting rid of http thing
                url = url[i+7:]
            elif url[i:i+8] == "https://":
                url = url[i+8:]
        for i in range(len(url)): # shaving off "www." and left
            if url[i:i+4] == "www.":
                url = url[i+4:]
                break
        for i in range(len(url)): # shaving off first '/' and right
            if url[i] == '/':
                url = url[:i]
                break
        for i in range(len(url) - 1, 0, -1): # backtracking to shave off TLD and right
            if url[i] == '.':
                url = url[:i]
                break
        for i in range(len(url) - 1, 0, -1): # removing rest of url from left
            if url[i] == '.':
                url = url[i+1:]
                break
        return url

    except IndexError: 
        print(f"Broken link; saving files to pid: /{os.getpid()}/")
        return os.getpid() # yes it's a cop out but at least it's unique

def main():
    # Get source code
    r = requests.get(sys.argv[1])

    # Handle to be able to parse source code
    soup = BeautifulSoup(r.text, 'html.parser')

    # All the hyperlink
    l_href = [link.get('href') for link in soup.find_all('a',href=True) ]

    # Use user-designated extension
    l_f90 = [link for link in filter(None,l_href) if link.endswith(sys.argv[3])]

    # Generate the full url
    final_url = websiteName(sys.argv[1])
    try: 
        os.mkdir(final_url)
        os.chdir(final_url)
    except:
        os.chdir(final_url)
        
    # Download everything and display
    l_url = [f"{sys.argv[2]}{name}" for name in l_f90]
    for name, url in zip(l_f90,l_url):
        name = fix(name)
        print (f'downloading {name} ({url})')
        r = requests.get(url)
        with open(f"{name}",'wb') as f:
            f.write(r.content)

    # Print new folder name
    print(f"done; results saved in /{final_url}/")

if __name__ == "__main__":
    main()
