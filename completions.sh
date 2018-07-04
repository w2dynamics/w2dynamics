#!/bin/bash
_hf_completionstr=${HF_PARAMETERS[@]}
_hf_longstr=${HF_LONGOPTS[@]}
_hf_longcomma=${_hf_longstr// /,}
_hf_shortstr=${HF_SHORTOPTS[@]}
_hf_shortjoined=${_hf_shortstr// /}
_hf_optionsstr="--${_hf_longstr// / --} -${_hf_shortstr// / -}"

_hf_complete() {
    local wordno=$COMP_CWORD
    local word=${COMP_WORDS[COMP_CWORD]}

    local opts=$(getopt -q -o "$_hf_shortjoined" -l "$_hf_longcomma" -- \
                 ${COMP_WORDS[@]})

    # Remove all options to get a consistent numbering for positional args
    for opt in $opts; do
        [ "$opt" == "--" ] && break
        let wordno--  # might then be less than zero if before positionals
    done

# Debugging
#    (
#       echo 
#	echo "${COMP_WORDS[@]:0:$COMP_CWORD+1}> ${COMP_WORDS[@]:$COMP_CWORD+1}"
#	echo "$wordno: $opts"
#    ) >>~/opt_debug

    if [ "$wordno" -le 1 -a "${word:0:1}" == "-" ]; then
        # Is an option because it is before a true quantity
        COMPREPLY=($(compgen -W "$_hf_optionsstr" -- "$word"))
    elif [ "$wordno" -eq 1 ]; then
        # Should be a HDF5 filename or a directory
        COMPREPLY=($(compgen -f -X "!*.hdf5" -- "$word")
                   $(compgen -d -X ".[^.]" -- "$word"))
    elif [ "$wordno" -eq 2 ]; then
        # Complete all the quantities which may be inside
        COMPREPLY=($(compgen -W "$_hf_completionstr" -- "$word"))
    fi
    return 0;
}
# FIXME: -o filenames works around a Debian bug causing appended spaces to 
# directories (different behaviour of complete and compgen). This poses problems
# in the unlikely case that there is a directory named like a quantity.
complete -o filenames -F _hf_complete hgrep

_dmft_completionstr=${DMFT_PARAMETERS[@]}
_dmft_complete() {
    local wordno=$COMP_CWORD
    local word=${COMP_WORDS[COMP_CWORD]}

    COMPREPLY=($(compgen -f -X "!*.in" -- "$word") 
               $(compgen -d -X ".[^.]" -S "/" -- "$word")
               $(compgen -W "$_dmft_completionstr" -S "=" -- "$word"))
    return 0;
}
complete -o nospace -F _dmft_complete DMFT.py
