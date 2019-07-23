# -*- coding: utf-8 -*-

from collections import namedtuple

Hit = namedtuple('Hit', 'Ref RefStart RefEnd Query QueryStart QueryEnd Score')
Seq = namedtuple('Seq', 'Name Seq Qual')
Segment = namedtuple('Segment', 'Start End Strand Len')
