import pandas as pd
class xls2dict:
    data = {}    
    def __init__(self, *source):
        
        def _prepare(wbk,msg):
            xls = pd.ExcelFile(wbk)
            for _ in xls.sheet_names:
                self.data[_] = pd.read_excel(xls, _)
            #print "%s into %s.data" % (msg,self.__class__.__name__)
            print msg
        if len(source):
            _prepare(source[0], "File '%s' readed successfully!" % source[0])
        else:
            from cStringIO import StringIO
            from IPython.display import display
            import fileupload 
            def _cbk(change):
                _prepare(StringIO(change['new']), "File '%s' uploaded successfully!" % change['owner'].filename)
            #
            _upload_widget = fileupload.FileUploadWidget()
            _upload_widget.observe(_cbk, names='data')
            display(_upload_widget)