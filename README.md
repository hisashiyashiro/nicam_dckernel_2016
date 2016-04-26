# nicam_dckernel_2016

* version: 0.1.0
* date: 2016/04/26

About NICAM dckernel package
--------------------------------
The nicam_dckernel_2016 is the kernel package of NICAM-DC.
NICAM-DC is the stand-alone package of the dynamical core of
Nonhydrostatic ICosahedral Atmospheric Model (NICAM).
The information of NICAM can be found at http://nicam.jp/



Installation
------------

#### Choose the installing platform.

In sysdep/ directory, there are Makedef.* files for different platforms.

If your installing platform is covered by one of these Makedef.* files,
then simply set NICAM_SYS environment variable using the chosen value,

If your testing platform is NOT covered by above Makedef.* files, please create the necessary Makedef files.
