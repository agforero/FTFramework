#!/usr/bin/env python3

import requests
import sys
from bs4 import BeautifulSoup

def fix(name):
    if len(name) == 0: return ""
    elif name[-1] == "/": return ""
    else: return str(fix(name[:-1])) + str(name[-1])

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
    l_url = [f"{sys.argv[2]}{name}" for name in l_f90]
    for name, url in zip(l_f90,l_url):
        name = fix(name)
        print (f'Downloading {name} ({url})')
        r = requests.get(url)
        with open(f"{name}",'wb') as f:
            f.write(r.content)

if __name__ == "__main__":
    main()