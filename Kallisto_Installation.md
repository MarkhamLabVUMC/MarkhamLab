# A Guide on how to Install Kallisto

Official Kallisto Download Commands: https://pachterlab.github.io/kallisto/download

**Windows**

The most straightforward way is to first download the Windows Subsystem for Linux (WSL). This lets developers install a Linux distribution (ex. Ubuntu) and use Linux applications, utilities, and Bash command-line tools directly on Windows.

1. Install WSL on your windows machine
    - Follow Microsoft's official guide (Link: https://learn.microsoft.com/en-us/windows/wsl/install)
2. Once WSL is set up and a Linux distribution is installed, open the Linux terminal
    - Type 'wsl' in the PowerShell/Windows Command Prompt and hit enter on your keyboard
3. For best practices, update the system to ensure you have the latest packages and security updates
    - Type 'sudo apt update' and 'sudo apt upgrade' one at a time in the PowerShell/Windows Command Prompt and hit enter on your keyboard
    - Note the 'sudo' command allows you to execute commands with superuser privileges
        - Might be prompted to enter your password
4. Download Kallisto from the GitHub Releases Page
    - Type 'wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz' in the PowerShell/Windows Command Prompt and hit enter on your keyboard
    - Note v0.46.1 is the latest version
5. Extract the Downloaded Archive
    - Type 'tar -xzvf kallisto_linux-v0.46.1.tar.gz' in the PowerShell/Windows Command Prompt and hit enter on your keyboard
6. Move Kallisto to your Path
    -After extraction, the Kallisto binary will be inside the extracted folder
    - For ease of use, move this binary to a location in your system's PATH
    - Common practice is to move it to '/usr/local/bin'
    - Type 'sudo mv kallisto_linux-v0.46.1/kallisto /usr/local/bin' in the PowerShell/Windows Command Prompt and hit enter on your keyboard
7. Verify Installation
    - Type 'kallisto' in the PowerShell/Windows Command Prompt and hit enter on your keyboard
    - Should display Kallisto help message, confirming it is installed and accessible from command line

**Mac**

The easiest way to download is via brew.

1. Install Homebrew on your mac machines
    - Follow Homebrew's official guide (Link: https://brew.sh/)
2. Install Kallisto
    - Type 'brew install kallisto' in the terminal and hit enter on your keyboard
3. Verify Installation
    - Type 'kallisto version' in the terminal and hit enter on your keyboard
    - Should display version number of Kallisto, confirming it is successfully installed
