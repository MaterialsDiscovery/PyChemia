from html.parser import HTMLParser
import sys


class MyHTMLParser(HTMLParser):
    """
    Create a subclass and override the handler methods
    """

    def error(self, message):
        pass

    def handle_starttag(self, tag, attrs):
        sys.stdout.write("'")

    def handle_endtag(self, tag):
        sys.stdout.write("'")

    def handle_data(self, data):
        sys.stdout.write(data)
