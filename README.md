# hifan--msutils
misc ms utilities, might get folded into a larger package at somepoint

# Installation
This package can be pip installed after you clone/download the repo. 
If you have git installed the recommended procedure is:

```
$ cd dir_for_repo_storage
$ git clone <repo_link_here>
$ cd hifan--mseutils
$ pip install -e .
$ python
>>> import mseutils
viola!
```

the `-e` flag on the pip install allows the code to be editible so if you pull new updates from this repository using `$ git pull` then the next time you import the code the  changes will be reflected. 


# Usage
I'll get around to describing usage soonish.. I'll make an example folder with some sample input data and a jupyter notebook with a few examples. 

# Notes
* This repo is a work in progress, there are no guarantees that the API won't change drastically, caveat emptor
* This isn't the speediest code in the world at the moment but it's designed to be easy to work with and play nicely with Jupyter Notebooks and such. 
    * hopefully that will lead to rapid development of future tools
    * mission critical sections can get re-worked in lower level numerical packages.