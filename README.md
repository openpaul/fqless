#fqless - FastQ Less
fqless is intendet to be a small less-like viewer for FastQ sequencing files. Allowing the user to get a quick glance at their sequecing data without the need of a heavy GUI/Java-bloatware.

It displays the Name and the color coded sequence. It hides the quality line from the user, as this is machine information for machines and not for human beeings.

It has no problem opening many gigabytes of fastq files, as it will not load everything, but will be buffering based on an index, it builds in memory. Thus after opening a large file, you will see hig CPU usage until the index is build, but you can already see everything.

fqless is released under the GPLv2 or any newer version of the GPL. It comes with no warranty. Use it as you wish.

## How to install

### On UNIX:
```
git clone ...
cd fqless
./configure
make
sudo make install
```

To uninstall just
```
sudo make uninstall
```

If you need or want to build your own Makefile run:
```
aclocal
autoconf
automake --add-missing
```



### On Windows:
Install UNIX then GOTO UNIX. 

Serioulsy: I dont have access to a Windows PC and no interest in porting this to windows. If you want to make a pull request with a working source.


## What it does
- It opens the file, shows the sequence color coded, so one can decide if the run quality is nice and fits the expectations.
It tries to guess the quality encoding of the file but you can also change that live.
- Allows you to convert files from one encoding into another
- Save FastQ as fasta (if someone wants that)

## What it does not
It does not allow to do anything else (no mapping, no trimming, no whatsoever).

KISS: Keep it simple, stupid.

## Features
- Parse FastQ [done]
- Detect encoding [todo]
- Color code sequence [done, could be better]
- allow switching between encodings [todo]
- allos conversion [todo]
- allow STDIN [yeah...]
- write fastq, allow conversion between encodings [todo]
- write fasta [todo]
