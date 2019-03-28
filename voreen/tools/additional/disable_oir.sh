#!/usr/bin/env bash
set -e

target=$1

geom_proc_props_path='/VoreenData/Workspace/ProcessorNetwork/Processors/Processor[@type="GeometryProcessor"]'

# Only modify workspaces that contain GeometryProcessor
if xmlstarlet sel -t -c $geom_proc_props_path $target > /dev/null; then

    if xmlstarlet sel -t -c $geom_proc_props_path'/Properties/Property[@mapKey="applyOrderIndependentTransparency_"]' $target > /dev/null; then
        xmlstarlet ed -P -L\
            --update $geom_proc_props_path'/Properties/Property[@mapKey="applyOrderIndependentTransparency_"]/@value' \
            --value "false"\
            $target
        echo "oir entry already present in $target, overwritten to false"
    else
        # Add an entry that disables oir
        xmlstarlet ed -P -L\
            --subnode $geom_proc_props_path'/Properties' -t elem -n Property \
            --subnode $geom_proc_props_path'/Properties/Property[not(@mapKey)]' -t attr -n mapKey -v "applyOrderIndependentTransparency_" \
            --subnode $geom_proc_props_path'/Properties/Property[not(@name)]' -t attr -n name -v "applyOrderIndependentTransparency_" \
            --subnode $geom_proc_props_path'/Properties/Property[not(@value)]' -t attr -n value -v "false" \
            $target

        echo "Added entry to disable oir in $target"
    fi
else
    echo "No GeometryProcessor present in $target"
fi
