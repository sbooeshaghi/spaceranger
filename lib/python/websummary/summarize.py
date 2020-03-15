"""Basic python interface for creating a websummary.

If run as "python summarize.py example" it
creates an example HTML websummary using the template in
example/summary.html and the data in example/data.json.
"""

import json
import os
import re
import sys

# Hard-coded template file to be used by `generate_html_summary`
DIST_TEMPLATE = 'dist/template.html'

FILES_TO_INLINE = {
    'tenx-websummary-script.min.js': 'dist/tenx-websummary-script.min.js',
    'tenx-websummary-styles.min.css': 'dist/tenx-websummary-styles.min.css',
}

# The regex to scan for includes into the summary template.  Syntax is
# [[ include LOCAL_FILE_PATH ]]
# File path can contain letters, digits, underscores, and dashes.  File path is local relative
# to the input template_dir.
TEMPLATE_REGEX = r"\[\[ include (?P<filename>[a-zA-Z.\/_\d-]+) \]\]"


def find_local_path(relative):
    """Compute path name of various resource files that are known relative to
    the location of this python file."""
    return os.path.join(os.path.dirname(__file__), relative)


def generate_html_summary(json_data, summary_contents, template_dir, output_filehandle, cls=None):
    """Generate an HTML summary.
    json_data:          a json-dumpable object containing data to be injected into the summary.
    summary_contents:   the inner HTML for the web summary
    template_dir:       a path to the folder containing the templates
    output_filehandle:  output where the final HTML will be written
    cls:                a JSONEncoder object to be passed to json.dumps(..., cls=cls)
    """
    template_file = DIST_TEMPLATE
    if template_dir and os.path.isfile(find_local_path(os.path.join(template_dir, 'template.html'))):
        template_file = find_local_path(os.path.join(template_dir, 'template.html'))
    with open(find_local_path(template_file), 'r') as template_fh:
        template = template_fh.read()

    for search_string, filename in FILES_TO_INLINE.items():
        with open(find_local_path(filename), 'r') as inline_fh:
            file_contents = inline_fh.read()
        template = template.replace(
            ('[[ %s ]]' % search_string), file_contents)

    # Replace inline subtemplates recursively from inner html.
    # We limit the number of loops to avoid infinite recursions.
    count = 0
    while True:
        if count > 100:
            raise ValueError("Maximum template recursion depth exceeded")
        count += 1
        match = re.search(TEMPLATE_REGEX, summary_contents)
        if match is None:
            break
        with open(os.path.join(template_dir, match.group("filename")), "r") as infile:
            file_contents = infile.read()
        summary_contents = summary_contents.replace(
            match.group(), file_contents)

    template = template.replace('[[ data.js ]]', json.dumps(json_data, cls=cls))
    template = template.replace('[[ summary.html ]]', summary_contents)

    output_filehandle.write(template)


def main(argv, outfile):
    """Executable entry point.

    Parses the command line, loads the data and content files, and writes
    the result to the given output.
    """
    json_file = "example/data.json"
    contents_file = "example/summary.html"
    if len(argv) > 1 and argv[1] == "cr":
        json_file = "tests/cr_tests/data/count_small.json"
        contents_file = "tests/cr_tests/summary.html"
    with open(find_local_path(json_file), 'r') as infile:
        data = json.load(infile)
    with open(find_local_path(contents_file), 'r') as infile:
        contents = infile.read()
    generate_html_summary(data, contents, os.path.join(
        os.path.dirname(__file__), 'example'), outfile)


if __name__ == "__main__":
    main(sys.argv, sys.stdout)
