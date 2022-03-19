from segments import fs


class NGSMutation(object):

    def __init__(self, mutType = None, site = None, size = None, substSize = None):
        self.mutType = mutType or fs.INVALID  # type FixedString
        self.site = site or 0  # type Integer
        self.size = size or 0  # type Integer
        self.substSize = substSize or 0  # type Integer

    @property
    def typeName(self):
        return "NGSMutation"

    @property
    def fsType(self):
        return fs.NGSMutation

    def defaultDict(self):
        return {
            'mutType' : self.mutType or fs.INVALID,
            'site' : self.site or 0,
            'size' : self.size or 0,
            'substSize' : self.substSize or 0,
        }

    def _description(self):
        return "NGSMutation: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return NGSMutation()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.mutType = json.get('mutType', fs.INVALID)
        self.site = json.get('site', 0)
        self.size = json.get('size', 0)
        self.substSize = json.get('substSize', 0)
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.mutType != None) if minimal else (self.mutType)): d['mutType'] = self.mutType
        if ((self.site != None) if minimal else (self.site)): d['site'] = self.site
        if ((self.size != None) if minimal else (self.size)): d['size'] = self.size
        if ((self.substSize != None) if minimal else (self.substSize)): d['substSize'] = self.substSize
        return d

class NGSRead(object):

    def __init__(self, sequenceId = None, success = None, overlap = None, mutations = None):
        self.sequenceId = sequenceId or 0  # type Integer
        self.success = success or False  # type Boolean
        self.overlap = overlap or 0  # type Integer
        self.mutations = mutations or []  # type [NGSMutation]

    @property
    def typeName(self):
        return "NGSRead"

    @property
    def fsType(self):
        return fs.NGSRead

    def defaultDict(self):
        return {
            'sequenceId' : self.sequenceId or 0,
            'success' : self.success or False,
            'overlap' : self.overlap or 0,
            'mutations' : self.mutations or [],
        }

    def _description(self):
        return "NGSRead: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return NGSRead()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.sequenceId = json.get('sequenceId', 0)
        self.success = json.get('success', False)
        self.overlap = json.get('overlap', 0)
        self.mutations = [ NGSMutation().loadFromJson(x, skipNull = skipNull) for x in json.get('mutations') or [] ]
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.sequenceId != None) if minimal else (self.sequenceId)): d['sequenceId'] = self.sequenceId
        if ((self.success != None) if minimal else (self.success)): d['success'] = self.success
        if ((self.overlap != None) if minimal else (self.overlap)): d['overlap'] = self.overlap
        if ((self.mutations != None) if minimal else (self.mutations)): d['mutations'] = [ x.json(skipTypes = skipTypes, minimal = minimal) for x in self.mutations ]
        return d

class NGSResult(object):

    def __init__(self, success = None, category = None, description = None, mutationCount = None, errorCount = None, targetSequence = None, targetDescription = None):
        self.success = success or False  # type Boolean
        self.category = category or ''  # type String
        self.description = description or ''  # type String
        self.mutationCount = mutationCount or 0  # type Integer
        self.errorCount = errorCount or 0  # type Integer
        self.targetSequence = targetSequence or ''  # type String
        self.targetDescription = targetDescription or ''  # type String

    @property
    def typeName(self):
        return "NGSResult"

    @property
    def fsType(self):
        return fs.NGSResult

    def defaultDict(self):
        return {
            'success' : self.success or False,
            'category' : self.category or '',
            'description' : self.description or '',
            'mutationCount' : self.mutationCount or 0,
            'errorCount' : self.errorCount or 0,
            'targetSequence' : self.targetSequence or '',
            'targetDescription' : self.targetDescription or '',
        }

    def _description(self):
        return "NGSResult: `{}`".format(", ".join([ "{}={}".format(k, v) for k, v in self.json(skipTypes = True).items() ]))

    def _newObjectOfSameType(self):
        return NGSResult()

    def clone(self):
        c = self._newObjectOfSameType()
        if hasattr(self, 'serialize'):
            c.deserialize(self.serialize())
        else:
            c.loadFromJson(self.json())
        return c

    def loadFromJson(self, json, skipNull = False):
        if not json:
            return self
        self.success = json.get('success', False)
        self.category = json.get('category', '')
        self.description = json.get('description', '')
        self.mutationCount = json.get('mutationCount', 0)
        self.errorCount = json.get('errorCount', 0)
        self.targetSequence = json.get('targetSequence', '')
        self.targetDescription = json.get('targetDescription', '')
        return self

    def json(self, skipTypes = False, minimal = False):
        d = { }
        if not skipTypes:
            d["type"] = self.typeName
        if ((self.success != None) if minimal else (self.success)): d['success'] = self.success
        if ((self.category != None) if minimal else (self.category)): d['category'] = self.category
        if ((self.description != None) if minimal else (self.description)): d['description'] = self.description
        if ((self.mutationCount != None) if minimal else (self.mutationCount)): d['mutationCount'] = self.mutationCount
        if ((self.errorCount != None) if minimal else (self.errorCount)): d['errorCount'] = self.errorCount
        if ((self.targetSequence != None) if minimal else (self.targetSequence)): d['targetSequence'] = self.targetSequence
        if ((self.targetDescription != None) if minimal else (self.targetDescription)): d['targetDescription'] = self.targetDescription
        return d
