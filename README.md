# MOLECULAR DYNAMICS SIMULATION

## Setup Instructions for Windows

To run the MOLECULAR DYNAMICS SIMULATION on a Windows system, you need to set up the following prerequisites:

1. **Cygwin**: Cygwin provides a large collection of GNU and Open Source tools that allow functionality similar to a Linux distribution on Windows. Follow the steps below to set it up:

   a. Download the Cygwin installer from the official website: https://www.cygwin.com/.
   
   b. Run the installer and select the appropriate options for your system. Make sure to include the required packages for C/C++ development.
   
   c. Complete the installation process by following the on-screen instructions.

2. **C Compiler**: You need a C compiler to build and execute the simulation. Cygwin provides the GCC compiler. Once you have installed Cygwin, the C compiler should be available automatically.

3. **Environment Variables**: To access the Cygwin and C environment from the Windows shell, you need to set the corresponding environment variables. Follow the steps below:

   a. Open the Control Panel on your Windows system.
   
   b. Go to "System and Security" and click on "System".
   
   c. Click on "Advanced system settings" on the left-hand side.
   
   d. In the System Properties window, click on the "Environment Variables" button.
   
   e. Under "User variables" or "System variables", find the "Path" variable and click on "Edit".
   
   f. Add the paths to your Cygwin installation and C compiler to the "Path" variable. For example:
      ```
      C:\cygwin64\bin\
      ```
      Make sure to adjust the paths based on your actual installation directories.
   
   g. Click "OK" to save the changes.

4. **PovRay**: PovRay is a ray-tracing software that may be required for rendering visualizations in the MOLECULAR DYNAMICS SIMULATION. Follow the steps below to install PovRay:

   a. Download the PovRay installer from the official website: https://www.povray.org/.
   
   b. Run the installer and follow the on-screen instructions to complete the installation process.
   
   c. Once installed, locate the PovRay installation directory.

5. **Environment Variable for PovRay**: To access PovRay from the Windows shell environment, you need to set the corresponding environment variable. Follow the steps below:

Add the Povray bin location to the environment variables path

'C:\Program Files\POV-Ray\v3.7\bin'

## Running the Simulation

Once you have completed the setup process, you are ready to run the MOLECULAR DYNAMICS SIMULATION on your Windows system. Make sure to open the Windows shell environment (e.g., Command Prompt, PowerShell) and navigate to the appropriate directory where the simulation code is located.
Run the Cygwin64 Terminal from the windows shell with the command

C:\cygwin64\bin\bash.exe

To compile and execute the simulation, use the following commands:

./build.sh
