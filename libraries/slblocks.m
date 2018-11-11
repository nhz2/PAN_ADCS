function blkStruct = slblocks
		% This function specifies that the library should appear
		% in the Library Browser
		% and be cached in the browser repository

		Browser.Library = 'collection_of_goodies';
		% 'collection_of_goodies' is the name of the library

		Browser.Name = 'PAN Library';
		% 'PAN Library' is the library name that appears 
             % in the Library Browser

		blkStruct.Browser = Browser; 