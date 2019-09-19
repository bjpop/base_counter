#!/bin/bash

set -e
errors=0

# Run unit tests
python base_counter/base_counter_test.py || {
    echo "'python python/base_counter/base_counter_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E base_counter/*.py || {
    echo 'pylint -E base_counter/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"
