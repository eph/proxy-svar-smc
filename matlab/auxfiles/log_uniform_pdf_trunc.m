function retval = log_uniform_pdf_trunc(x, scale)
retval = 0;
if (x > scale) 
    retval= -10000000;
end

end 