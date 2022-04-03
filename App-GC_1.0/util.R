choose_directory = function(caption = 'Select data directory') {
  if (exists('choose.dir')) {
    utils::choose.dir('', caption = caption) 
  } else {
    tcltk::tk_choose.dir('', caption = caption)
  }
}