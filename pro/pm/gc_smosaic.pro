pro gc_smosaic
name= ["NGC 104", $
       "NGC 288", $
       "NGC 362", $
       "NGC 1851", $
       "NGC 1904", $
       "NGC 2298", $
       "NGC 4147", $
       "NGC 4590", $
       "NGC 5024", $
       "NGC 5139", $
       "NGC 5272", $
       "NGC 5466", $
       "NGC 5897", $
       "NGC 5904", $
       "NGC 6093", $
       "NGC 6121", $
       "NGC 6144", $
       "NGC 6171", $
       "NGC 6205", $
       "NGC 6218", $
       "NGC 6254", $
       "NGC 6341", $
       "NGC 6362", $
       "NGC 6397", $
       "NGC 6584", $
       "NGC 6626", $
       "NGC 6656", $
       "NGC 6712", $
       "NGC 6752", $
       "NGC 6779", $
       "NGC 6809", $
       "NGC 6838", $
       "NGC 7078", $
       "NGC 7089", $
       "NGC 7099", $
       "Pal 3", $
       "Pal 5"]
for ii=0,n_elements(name)-1 do begin
    hogg_ned_smosaic, name[ii],1024,768,nskyframe=5
endfor
return
end
