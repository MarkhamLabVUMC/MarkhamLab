# External Hard Drive Guidelines

This provides key points to implement when working with the drive (Windows) to ensure everything goes smooth.

**Assigning a Letter**

Sometimes, a letter is not assigned to the external hard drive, so it won't show up in the File Explorer. You need to manually assign a letter by following these steps:

1. Right-click on the start button and select "Disk Management"
2. Find the external hard drive in the list
3. Right-click on the drive and choose "Change Drive Letter and Paths"
4. Assign a drive letter that is not in use (Highly recommend to use 'E' as all code assumes that is the letter assigned to the drive)

**Setup/Mounting**

Since the data is on external hard drive, you have to mount it in wsl to be able to use it. Follow these steps:

1. Open WSL
    - Open Windows Terminal or Powershell and type 'wsl'
2. Create a Mount Point
    - Decide where you want to mount your network drive in your WSL file system
    - Common choice is under /mnt/, like '/mnt/e'
    - Create mount point by typing 'mkdir -p /mnt/e' in the terminal and press enter
3. Mount the Network Drive
    - Type 'sudo mount -t drvfs '\\Book-m08q8a5m4r\e' /mnt/e' in the terminal and press enter.
4. Verify the Mount
    - Type 'ls /mnt/e' in the terminal and press enter
    - If you see the contents then it was successful 

*Note this is a one time item to do for your machine.

**Ejecting**

Ejecting the hard drive safely is crucial so the files do not go corrupt. You cannot simply physically remove the USB connection until you eject it by clicking on the USB icon on the bottom right corner. However, this method doesn't always work. So, here is what you can do:

1. Open Command Prompt as an administrator
    - Right-click on the app and press 'Run as Administrator'
2. Type 'diskpart' and press enter
3. Type 'list volume' and press enter
4. Find the volume number of your external hard drive
5. Type 'select volume [number]', where '[number]' is the number of the drive and press enter
6. Type 'remove all dismount' and press enter
7. Then, you can physically disconnect the external hard drive from the machine