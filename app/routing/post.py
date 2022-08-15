"""
Created by Daniel Fu (Reif Lab, Duke University) at 4/10/2021

Name        : post
Project     : cadaxisdna
Description : post-processing heuristics
Interpreter : Python 3.7.4
"""

from .helper import log, strandnav


class Mixin:
    def rebalance_short_strands(self):
        """
        Check for short strands
        :return:
        """
        for staple in self.get_all_staples():
            if len(staple) < 15:
                log.system("Found short staple ({}): {}".format(len(staple), [n.numid for n in staple]))
                log.system("Attempting to rebalance nicks with its adjacent strands.")
                strand5 = strandnav.get5strand(staple)
                strand3 = strandnav.get3strand(staple)
                log.system("5' end ({}): {}".format(len(strand5), [n.numid for n in strand5]))
                log.system("3' end ({}): {}".format(len(strand3), [n.numid for n in strand3]))
                self.balance(staple)


if __name__ == "__main__":
    pass
