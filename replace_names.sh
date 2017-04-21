#!/usr/bin/env bash

while read from to; do
  sed -ie "s/$from/$to/" pop_names.csv ;
done < truncated.list