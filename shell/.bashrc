# Sets up BASH with the CoDeS3D environment.
# First, source the user's default environment
test -e "/etc/bash.bashrc" && source "/etc/bash.bashrc"
test -e "~/.bashrc" && source "~/.bashrc"
test -e "~/.bash_profile" && source "~/.bash_profile"

# Then, add the CoDeS3D prompt
PS1="$PS1 CoDeS3D > "
