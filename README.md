Build with cargo build. This will also create a bash completion script. Something like:

```
warning: completion file is generated: "/home/niklas/code/rust/target/debug/build/seq_tools-85cd6d057134563e/out/seq_tools.bash"
```

This file contains a script to make tab-completion work. To install the script, copy it to:

```
~/.local/share/bash-completion/completions/seq_tools.bash
```

If that directory does not exist, you can create it yourself.
