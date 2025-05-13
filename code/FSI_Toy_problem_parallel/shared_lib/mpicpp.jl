# Load the module and generate the functions
module mpicpp
  using CxxWrap
  @wrapmodule(() -> joinpath("/path/to/shared_lib","libmpicpp"))

  function __init__()
    @initcxx
  end
end
