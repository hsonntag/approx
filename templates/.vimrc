 set nocompatible
 
 set autoindent
 set smartindent

 set tabstop=4
 set shiftwidth=4

 set showmatch

 set vb t_vb=

 set ruler

 set nohls

 set incsearch

 " set virtualedit=all

 filetype plugin on
 set grepprg=grep\ -nH\ $*
 filetype indent on
 let g:tex_flavor='latex'

 colorscheme pablo
 syntax on
 set clipboard+=unnamed
 set mouse=a
 set list " we do what to show tabs, to ensure we get them
          " out of my files
 set listchars=tab:>-,trail:- " show tabs and trailing
 set expandtab " no real tabs please!
 set foldenable " Turn on folding
 set foldmarker={,}
