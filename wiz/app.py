from flask import Flask, render_template, request

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def form():
    return render_template('angular.html')

# @app.route('/hello', methods=['GET', 'POST'])
# def hello():
#     return render_template('greeting.html', data=request.form)

# @app.route('/<command>', methods=['GET'])
# def render_command_opts(command):
#     return render_template(command)

def main():
    app.run(debug=True, host='0.0.0.0')

if __name__ == "__main__":
    main()
