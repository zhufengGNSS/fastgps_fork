Below is a brief summary of how to build the fastgps receiver under Windows.
I've piled all the new files in the same directory as the trunk folder of the fastgps project.
Which is where the song and dance below should be executed from. 

The following two bakefiles have been included for building the fastgps software receiver under MicroSoft Windows (MSW).

fastgps_windows.bkl
wxfastgps_windows.bkl

The first builds a console application, while the second the simple wxWidgets version (which calls wx.bkl, wx_win32.bkl and wx_unix.bkl presets).

These two files can be used to generate native makefiles under MSW using the free bakefile program.

http://www.bakefile.org/index.html

The actual makefiles are generated with the following lines,

bakefile -f mingw -I.. fastgps_windows.bkl
bakefile -f mingw -I.. wxfastgps_windows.bkl

or

bakefile -f borland -I.. fastgps_windows.bkl
bakefile -f borland -I.. wxfastgps_windows.bkl

both of these lines will generate/overwrite the file makefile.gcc or makefile.bcc
For other compilers, substitude the name of the compiler above, ex. Watcom is watcom, MS VC++ is msvc.

This file can then be used to build an executable with the following line,
make -f filename  

(be sure you use the make command specific for each compiler, i.e. don't run the mingw make on the borland makefile!)

Pre generated makefiles for the console and wxWidgets versions have been included; fastgps_makefile.gcc and wxfastgps_makefile.gcc.

The tools you will need are:
    1) A C++ compiler, accessable from the command line (i.e. in the system PATH). The above examples were tested using the mingw
    gcc compiler (http://www.mingw.org/) and the Borland free C++ compiler (http://www.codegear.com/downloads/free/cppbuilder)

    2) A make utility, often included with compiler, callable from the MS-DOS command line. There are several of these available
    free on the internet, including; MSYS shell environment which comes with MingW (http://www.mingw.org/), cygwin and GNU make. 

    3) For the wxWidgets build, you will need the wxWidgets libraries installed and in your system PATH. Please see 
    http://www.wxwidgets.org/ or http://wxpack.sourceforge.net/ for details on installing wxWidgets libraries. wxWidgets is 
    compiler dependant, so for the two examples above you would need both the gcc and bcc libraries built. Check the wxWidgetsX.X\lib directory

NOTE: If you don't care to use makefiles, its always possible to copy all of the project files into a separate directory 
and use any number of C++ IDE's (such as Code::Blocks, http://www.codeblocks.org/) to create a new project, which would make 
make available the debug tools contained in the IDE. 

a few things to check if it doesn't work,

1) check your compiler is setup correctly. For example, at the command line type "gcc" or "bcc32". If it complains about 
no files then it is OK, if it can't find the command then there is a compiler setup problem.
2) Paths, paths, paths. Make sure the compiler can see the libraries and header files. If it is moaning about not being 
able to find a standard header file like stdio.h, the compiler paths are not right. 
3) If it complains about not being able to find stdint.h, this means you may need to change the types.h file to include your
compiler predefine (see http://predef.sourceforge.net/precomp.html). It is assuming the unix data types are available,
which they are for the mingw compiler but not for most others.

I realize it is difficult to generalize makefiles on Windows, so I'm sure there are conditions that will derail this plan.
So if you are able, I'd consider switching to Linux (or MacOS) where portable makefiles are alot easier :-)
If anyone wants to write up a better guide, be my guest!

Cheers, STG

just in case people find this mess useful,
License: GNU GPL 
