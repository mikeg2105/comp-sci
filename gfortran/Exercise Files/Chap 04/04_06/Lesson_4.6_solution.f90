module todo_m
    implicit none
    private
    public :: interact, read_previous_tasks, MAX_TASKS, TASK_LENGTH

    integer, parameter :: MAX_TASKS = 99
    integer, parameter :: TASK_LENGTH = 100
    character(len=100) :: todo_file
contains
    subroutine add(tasks, num_tasks)
        character(len=TASK_LENGTH), intent(inout) :: tasks(MAX_TASKS)
        integer, intent(inout) :: num_tasks

        if (num_tasks == MAX_TASKS) then
            print *, "Can't save any more tasks."
        else
            print *, "What's the task?"
            num_tasks = num_tasks + 1
            read(*, '(A)') tasks(num_tasks)
        end if
    end subroutine

    subroutine ask_user(response)
        character(len=1), intent(out) :: response

        print *, "What would you like to do? (a)dd, (d)elete, (q)uit"
        read(*, '(A)') response
    end subroutine

    subroutine delete(tasks, num_tasks)
        character(len=TASK_LENGTH), intent(inout) :: tasks(MAX_TASKS)
        integer, intent(inout) :: num_tasks

        integer :: i
        integer :: status
        integer :: task

        if (num_tasks > 0) then
            do
                print *, "Which one would you like to delete?"
                read(*, *, iostat = status) task
                if (status == 0) then
                    if (task < 1) then
                        print *, "Task number must be >= 1."
                    else if (task > num_tasks) then
                        print *, "Task number must be <= ", num_tasks
                    else
                        do i = task, num_tasks-1
                            tasks(i) = tasks(i+1)
                        end do
                        num_tasks = num_tasks - 1
                        exit
                    end if
                else
                    print *, "Sorry, I didn't understand that."
                end if
            end do
        else
            print *, "No tasks to delete."
        end if
    end subroutine

    subroutine determine_file_name
        select case (command_argument_count())
        case (0)
            todo_file = "todo.txt"
        case (1)
            call get_command_argument(1, todo_file)
        case default
            print *, "Too many command line arguments."
            print *, "Usage: todo [file]"
            print *, "  where file is optional, and is 'todo.txt' by default"
            stop 1
        end select
    end subroutine

    subroutine interact(tasks, num_tasks)
        character(len=TASK_LENGTH), intent(inout) :: tasks(MAX_TASKS)
        integer, intent(inout) :: num_tasks

        character(len=1) :: response

        do
            call report_tasks(tasks, num_tasks)
            call ask_user(response)
            select case (response)
            case ('a')
                call add(tasks, num_tasks)
            case ('d')
                call delete(tasks, num_tasks)
            case ('q')
                exit
            case default
                print *, "Sorry, I didn't understand that."
            end select
            call save(tasks, num_tasks)
        end do
    end subroutine

    subroutine read_previous_tasks(tasks, num_tasks)
        character(len=TASK_LENGTH), intent(out) :: tasks(MAX_TASKS)
        integer, intent(out) :: num_tasks

        call determine_file_name
        call read_tasks_from_file(tasks, num_tasks)
    end subroutine

    subroutine read_tasks_from_file(tasks, num_tasks)
        character(len=TASK_LENGTH), intent(out) :: tasks(MAX_TASKS)
        integer, intent(out) :: num_tasks

        integer :: i
        integer :: status
        integer :: unit

        open(newunit=unit, file=trim(todo_file))
        do i = 1, MAX_TASKS
            read(unit, '(A)', iostat=status) tasks(i)
            if (status /= 0) then
                num_tasks = i - 1
                exit
            else
                num_tasks = i
            end if
        end do
        close(unit)
    end subroutine

    subroutine report_tasks(tasks, num_tasks)
        character(len=TASK_LENGTH), intent(in) :: tasks(MAX_TASKS)
        integer, intent(in) :: num_tasks

        integer :: i

        if (num_tasks > 0) then
            print *, "Here are your current tasks."
            do i = 1, num_tasks
                print '(I3,") ",A)', i, trim(tasks(i))
            end do
        else
            print *, "You have no tasks."
        end if
    end subroutine

    subroutine save(tasks, num_tasks)
        character(len=TASK_LENGTH), intent(in) :: tasks(MAX_TASKS)
        integer, intent(in) :: num_tasks

        integer :: i
        integer :: unit

        open(newunit=unit, file=trim(todo_file), status="replace")
        do i = 1, num_tasks
            write(unit, '(A)') trim(tasks(i))
        end do
        close(unit)
    end subroutine
end module

program todo
    use todo_m, only: interact, read_previous_tasks, MAX_TASKS, TASK_LENGTH
    implicit none

    integer :: num_tasks
    character(len=TASK_LENGTH) :: tasks(MAX_TASKS)

    call read_previous_tasks(tasks, num_tasks)
    call interact(tasks, num_tasks)
end program
