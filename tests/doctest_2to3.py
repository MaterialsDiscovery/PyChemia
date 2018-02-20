import re
import sys
import doctest


class DocTestChecker(doctest.OutputChecker):
    def check_output(self, want, got, optionflags):
        if sys.version_info[0] > 2:
            want = re.sub("u'(.*?)'", "'\\1'", want)
            want = re.sub('u"(.*?)"', '"\\1"', want)
        else:
            want = re.sub("b'(.*?)'", "'\\1'", want)
            want = re.sub('b"(.*?)"', '"\\1"', want)
        return doctest.OutputChecker.check_output(self, want, got, optionflags)


def doctest_suite(mod):
    return doctest.DocTestSuite(mod, checker=DocTestChecker())
