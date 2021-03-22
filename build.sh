pyinstaller probeSearch.py -F \
    --name "ProbeSearch" \
    --icon="icons/myicon.ico" \
    --add-binary='/System/Library/Frameworks/Tk.framework/Tk':'tk' \
    --add-binary='/System/Library/Frameworks/Tcl.framework/Tcl':'tcl' \
    --add-data="icons/*.svg:icons/." \
    --add-data="*.ui:." \
    --add-data="*.yml:." \
    --clean