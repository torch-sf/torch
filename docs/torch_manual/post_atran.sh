#!/usr/bin/env bash

# Meant for use by A Tran only, to post to web
# in order to prevent search engine crawlers it is weakly password protected
rsync -avh main.pdf pluto:~/www/torch/manual.pdf
