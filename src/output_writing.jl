#file containing outputting functions


function write_output(wavi::AbstractModel,output::Output,clock::Clock)
    #loop over every entry in the dictionary, for each entry, get the corresponding matrix
    for (key,val) in output.out_dict
        field = fetch_val(wavi,val)
        save_field(output,clock)
    end
    return nothing
end


function fetch_val(wavi,val::Expr)


    return nothing
end