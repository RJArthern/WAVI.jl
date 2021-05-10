#file containing outputting functions


function write_output(wavi::AbstractModel)
    #initilize empty dictionary for output
    dict = Dict()
    #loop over every entry in the dictionary, for each entry, get the corresponding matrix
    for (key,val) in wavi.output.out_dict
        field = fetch_val(wavi,val)
        dict[key] = field
    end
    return dict
end


function fetch_val(wavi_,val::Expr)
    #define a function to eval
    #wavi = wavi
    println(val)
    
    field_fn = @eval f(wavi)=$val
            
    field = field_fn()

    return field
end

function initialize_output_variables!(wavi)
    #for each value, evaluate that expression now
    dict = Dict()
    
    for (key,val) in wavi.output.out_dict
        println(eval(val))
    end
    wavi = @set wavi.output.out_dict_values = dict
    
    return nothing
end